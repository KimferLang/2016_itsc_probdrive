module FeaturesetExtractor

using DataFrames
using HDF5, JLD

using ProbDrive2016ITSC.Curves
using ProbDrive2016ITSC.Trajdata
using ProbDrive2016ITSC.StreetNetworks
using ProbDrive2016ITSC.Features
include("../io/filesystem_utils.jl")

export
    CSVFileSet,

    FEATURES,
    BEHAVIORS,

    get_validfind_regions,
    merge_region_segments,

    gen_featureset_from_regions,
    gen_featureset_from_validfinds,
    extract_csvfile_set,
    extract_csvfile_sets,

    does_violate_filter,

    create_dataframe_with_feature_columns


type CSVFileSet
    carid::Int
    csvfile :: AbstractString
    streetmapbasename :: AbstractString

    only_extract_these_frame_ranges::Vector{Int}

    lanechanges_normal    :: Vector{Int}
    lanechanges_postpass  :: Vector{Int}
    lanechanges_arbitrary :: Vector{Int}
    carfollow             :: Vector{Int}
    freeflow              :: Vector{Int}
end

const FEATURES = allfeatures()
const BEHAVIORS = ["lanechange", "carfollow", "freeflow"]

function extract_all_lanechanges(
    pdset::PrimaryDataset,
    sn::StreetNetwork,
    carid::Integer,
    validfind_contains_carid::AbstractVector{Bool},
    )

    #=
    Returns a list of frames in which the vehicle is, for the first time, in a new lane
    =#

    lanechange_validfinds = Int[]


    for validfind = 2 : length(validfind_contains_carid)
        if validfind_contains_carid[validfind] && validfind_contains_carid[validfind-1]

            cur_validfind = validfind-1
            cur_frameind = validfind2frameind(pdset, cur_validfind)
            cur_carind = carid2ind(pdset, carid, cur_validfind)

            fut_validfind = validfind
            fut_frameind = validfind2frameind(pdset, fut_validfind)
            fut_carind = carid2ind(pdset, carid, fut_validfind)

            # if isna(get(pdset, :lanetag, carind, cur_frameind, cur_validfind))
            #     println("validfind_contains_carid: ", validfind_contains_carid[validfind])
            #     println("validfind_contains_carid: ", validfind_contains_carid[cur_validfind])
            #     println("idinframe: ", idinframe(pdset, carid, validfind))
            #     println("idinframe: ", idinframe(pdset, carid, cur_validfind))
            #     matind = pdset.dict_other_idmap[carid]
            #     println("matind: ", matind)
            #     println("mat_other_indmap: ", pdset.mat_other_indmap[validfind, matind])
            #     println("mat_other_indmap: ", pdset.mat_other_indmap[cur_validfind, matind])
            #     println("validfind: ", validfind)
            #     println("cur_validfind: ", cur_validfind)
            #     println("cur_frameind: ", cur_frameind)
            #     println("carind: ", carind)
            #     println(get(pdset, :lanetag, carind, cur_frameind, cur_validfind))
            #     println("\n")
            #     # println(pdset.df_other[cur_validfind, :])
            # end

            fut_lanetag = get(pdset, :lanetag, fut_carind, fut_frameind, fut_validfind)::LaneTag
            cur_lanetag = get(pdset, :lanetag, cur_carind, cur_frameind, cur_validfind)::LaneTag
            cur_lane = get_lane(sn, cur_lanetag)



            if cur_lanetag != fut_lanetag
                # potential lane change - need to verify it wasn't just us crossing a tile section
                # or road segment

                lanechange = same_tile(cur_lanetag, fut_lanetag) || !has_next_lane(sn, cur_lane)
                if !lanechange
                    cur_lane = next_lane(sn, cur_lane)
                    cur_lanetag = cur_lane.id
                    lanechange = fut_lanetag != cur_lanetag
                end

                if lanechange
                    push!(lanechange_validfinds, fut_validfind)
                end
            end
        end
    end

    lanechange_validfinds
end
function identify_lane_changes!(
    validfind_is_lanechange::AbstractVector{Bool},
    validfind_contains_carid::AbstractVector{Bool},
    lanechange_validfinds::AbstractVector{Int},
    pdset::PrimaryDataset,
    carid::Integer,
    lat_vel_threshold::Real=0.1, # [m/s]
    )

    #=
    scans back and fore for each lane change in lanechange_validfinds,
    marking all validfinds with an absolute lateral velocity greater than the threshold
    =#

    @assert(length(validfind_is_lanechange) == length(validfind_contains_carid))

    fill!(validfind_is_lanechange, false)

    for lanechange_validfind in lanechange_validfinds
        @assert(validfind_contains_carid[lanechange_validfind])
        validfind_is_lanechange[lanechange_validfind] = true

        # scan forward
        validfind = lanechange_validfind
        while validfind_inbounds(pdset, validfind+1) &&
              validfind_contains_carid[validfind+1] &&
              abs(get(pdset, :velFy, carid2ind(pdset, carid, validfind+1), validfind+1)) > lat_vel_threshold

            validfind += 1
            validfind_is_lanechange[validfind] = true
        end

        # scan backward
        validfind = lanechange_validfind
        while validfind_inbounds(pdset, validfind-1) &&
              validfind_contains_carid[validfind-1] &&
              abs(get(pdset, :velFy, carid2ind(pdset, carid, validfind-1), validfind-1)) > lat_vel_threshold

            validfind -= 1
            validfind_is_lanechange[validfind] = true
        end
    end

    validfind_is_lanechange
end

function extract_csvfile_sets(
    csvfile::AbstractString,
    streetmapbasename::AbstractString,
    pdset::PrimaryDataset,
    sn::StreetNetwork,
    )


    #=
    Loads the given pdset and attempts to automatically classify segments into:
     - lane change:  absolute lateral velocity is > 0.1 m/s and lane centerline switch occurs
        - normal:    lane change using car follow
        - postpass:  used to be in the target lane behind a car and am now passing in front of said car
        - arbitrary: not normal or postpass
     - free flow:    front timegap > 3s or front vehicle faster by 0.5 m/s
     - car follow:   not in free flow

     NOTE(tim): for now it classifies everything as a 'normal' lanechange
    =#

    retval = CSVFileSet[]

    N = nvalidfinds(pdset)
    validfind_contains_carid = falses(N)
    validfind_is_lanechange = falses(N)
    validfind_is_freeflow = falses(N)

    for carid in get_carids(pdset)
        push!(retval, extract_csvfile_set(carid, csvfile, streetmapbasename, pdset, sn,
                validfind_contains_carid, validfind_is_lanechange, validfind_is_freeflow
            ))
    end
    retval
end
function extract_csvfile_sets(
    csvfile::AbstractString,
    streetmapbasename::AbstractString,
    )

    pdset = load(joinpath(PRIMARYDATA_DIR, toext("primarydata_" * csvfile, "jld")))["pdset"]
    sn = load(STREETMAP_BASE*streetmapbasename*".jld")["streetmap"]
    extract_csvfile_sets(csvfile, streetmapbasename, pdset, sn)
end

function get_validfind_regions(behavior::AbstractString, csvfileset::CSVFileSet)

    # in(behavior, BEHAVIORS) || error("unknown behavior $behavior")

    if behavior == "lanechange"
        merge_region_segments([csvfileset.lanechanges_normal, csvfileset.lanechanges_postpass, csvfileset.lanechanges_arbitrary])
    elseif behavior == "carfollow"
        merge_region_segments([csvfileset.carfollow])
    elseif behavior == "freeflow"
        merge_region_segments([csvfileset.freeflow])
    else behavior == "all"
        merge_region_segments([csvfileset.lanechanges_normal, csvfileset.lanechanges_postpass,
                               csvfileset.lanechanges_arbitrary, csvfileset.freeflow,
                               csvfileset.carfollow])
    end
end

function is_strictly_monotonically_increasing(vec::AbstractVector{Int})
    if isempty(vec)
        return true
    end

    val = vec[1]
    for i = 2 : length(vec)
        val2 = vec[i]
        if val2 ≤ val
            return false
        end
        val = val2
    end
    true
end
function is_strictly_monotonically_nondecreasing(vec::AbstractVector{Int})
    if isempty(vec)
        return true
    end

    val = vec[1]
    for i = 2 : length(vec)
        val2 = vec[i]
        if val2 < val
            return false
        end
        val = val2
    end
    true
end
function merge_region_segments(vec::AbstractVector{Int})

    #=
    sorts a set of region segments to be in monotonically increasing order
    will also merge them as necessary if they overlap
    =#

    if isempty(vec)
        return Int[]
    end

    n = length(vec)
    @assert(mod(n,2) == 0)

    validfind_min, validfind_max = extrema(vec)

    validfinds = falses(validfind_max-validfind_min+1)
    for i = 1 : div(n,2)
        validfind_lo = vec[2i-1]-validfind_min+1
        validfind_hi = vec[2i]-validfind_min+1
        for j = validfind_lo : validfind_hi
            validfinds[j] = true
        end
    end

    retval = Int[]
    validfind_lo = findfirst(validfinds)
    while validfind_lo != 0
        validfind_hi = findnext(val->!val, validfinds, validfind_lo+1)
        if validfind_hi == 0
            validfind_hi = length(validfinds)
        end
        push!(retval, validfind_lo+validfind_min-1)
        push!(retval, validfind_hi+validfind_min-1)

        validfind_lo = findnext(validfinds, validfind_hi+1)
    end

    retval
end
function calc_row_count_from_region_segments(validfind_regions::AbstractVector{Int})

    #=
    Computes the number of frames encoded in validfind_regions

    validfind_regions is an array in which subsequent value pairs correspond to continuous frames
    containing a particular behavior.
    validfind_regions should be strictly monotonically increasing

    ex: [1,4,8,10] means 1→4 and 8→10 are continuous frames
    =#

    n = length(validfind_regions)
    if n == 0
        return 0
    end

    if !is_strictly_monotonically_nondecreasing(validfind_regions)
        println("not strincy increasing: ", validfind_regions)
    end

    @assert(mod(n, 2) == 0) # must be an even number
    @assert(is_strictly_monotonically_nondecreasing(validfind_regions))

    estimated_row_count = 0
    index_lo = 1
    while index_lo < length(validfind_regions)
        estimated_row_count += validfind_regions[index_lo+1] - validfind_regions[index_lo] + 1
        index_lo += 2
    end
    estimated_row_count
end
function create_dataframe_with_feature_columns{F<:AbstractFeature}(features::AbstractVector{F}, nrows::Int)
    df = DataFrame()
    df[:validfind] = DataArray(Int, nrows)
    for f in features
        df[symbol(f)] = DataArray(Float64, nrows)
    end
    df
end

end # module