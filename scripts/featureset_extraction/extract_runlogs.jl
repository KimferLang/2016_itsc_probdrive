using ProbDrive2016ITSC
using ProbDrive2016ITSC.StreetNetworks.RoadNetwork

include(Pkg.dir("ProbDrive2016ITSC", "src", "io", "filesystem_utils.jl"))

const STREETNET_CACHE = Dict{AbstractString, StreetNetwork}()
const PRIMARYDATA_DIR = "/media/tim/DATAPART1/Data/Bosch/processed/primarydata/"
const OUTPUT_DIR = "/media/tim/DATAPART1/PublicationData/2015_TrafficEvolutionModels/realworld/"
const RUNLOG_OUTPUT_DIR = joinpath(OUTPUT_DIR, "runlogs")
const RUNLOG_DIR = joinpath(OUTPUT_DIR, "runlogs")


const CSVFILESETS = (
            CSVFileSet(RunLogs.ID_EGO,
                        "path_to_runlog.csv",
                        "???",
                        Int[],
                        Int[], # lanechange normal
                        Int[], # lanechanges postpass
                        Int[], # lanechanges arbitrary
                        Int[], # car follow
                        Int[], # free flow
                        ),
            ... # MOAR HERE
        )

#############################################

@assert(isdir(OUTPUT_DIR))
if !isdir(RUNLOG_OUTPUT_DIR)
    mkdir(RUNLOG_OUTPUT_DIR)
end

#############################################

extract_params = PrimaryDataExtractionParams()

tot_frames = 0
tot_freeflow = 0
tot_following = 0
tot_lanechange = 0

tic()
for (csvfileset_index, csvfileset) in enumerate(CSVFILESETS)

    println("csvfileset ", csvfileset_index, " / ", length(CSVFILESETS))
    csvfilename = csvfileset.csvfile
    csvfilebase = basename(csvfilename)
    csvfilebase_noext = splitext(csvfilebase)[1]

    header, trajdata, sn = load_header_trajdata_and_streetmap(csvfilename)
    # extract_params.csvfileset = csvfileset
    extract_params.frameinds = csvfileset.only_extract_these_frame_ranges
    runlogs = extract_runlogs(trajdata, sn, extract_params, header)::AbstractVector{RunLog}

    for (i,runlog) in enumerate(runlogs)
        runlogname = joinpath(RUNLOG_OUTPUT_DIR, @sprintf("primarydata_%s_%d.jld", splitext(csvfilebase)[1], i))
        JLD.save(runlogname, "runlog", runlog)
    end

    for runlog in runlogs

        tot_frames += nframes(runlog)
        for frame in 1 : nframes(runlog)
            colset = colset2id(runlog, ID_EGO, frame)
            tot_freeflow += is_behavior_flag_set(runlog, colset, frame, ContextClass.FREEFLOW)
            tot_following += is_behavior_flag_set(runlog, colset, frame, ContextClass.FOLLOWING)
            tot_lanechange += is_behavior_flag_set(runlog, colset, frame, ContextClass.LANECHANGE)
        end
    end
end
toc()

println("tot_frames:     ", tot_frames)
println("tot_freeflow:   ", tot_freeflow)
println("tot_following:  ", tot_following)
println("tot_lanechange: ", tot_lanechange)