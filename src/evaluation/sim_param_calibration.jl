export
        ParamsHistobin,

        calc_trace_deviation,

        allocate_empty_histobin,
        update_histobin!,
        calc_histobin,

        KL_divergence_categorical,
        KL_divergence_dirichlet,
        KL_divergence_univariate_gaussian

type ParamsHistobin
    discx :: LinearDiscretizer
    discy :: LinearDiscretizer

    ParamsHistobin(discx::LinearDiscretizer, discy::LinearDiscretizer) = new(discx, discy)
    function ParamsHistobin(binedgesx::Vector{Float64}, binedgesy::Vector{Float64})
        discx = LinearDiscretizer(binedgesx)
        discy = LinearDiscretizer(binedgesy)
        new(discx, discy)
    end
end

function calc_trace_deviation(pdset::PrimaryDataset, sn::StreetNetwork, seg::PdsetSegment)

    initial_carind = carid2ind(pdset, seg.carid, seg.validfind_start)
    velFx = get(pdset, :velFx, initial_carind, seg.validfind_start)
    velFy = get(pdset, :velFy, initial_carind, seg.validfind_start)
    initial_speed = hypot(velFx, velFy)

    final_carind = carid2ind(pdset, seg.carid, seg.validfind_end)
    posGx_A = get(pdset, :posGx, initial_carind, seg.validfind_start)
    posGy_A = get(pdset, :posGy, initial_carind, seg.validfind_start)
    posGx_B = get(pdset, :posGx, final_carind, seg.validfind_end)
    posGy_B = get(pdset, :posGy, final_carind, seg.validfind_end)

    Δs, Δd = frenet_distance_between_points(sn, posGx_A, posGy_A, posGx_B, posGy_B)
    if isnan(Δs)
        warn("NaN Δs: ", (Δs, Δd))
        println((posGx_A, posGy_A, posGx_B, posGy_B))
        println(project_point_to_streetmap(posGx_A, posGy_A, sn))
    end
    @assert(!isnan(Δd))
    @assert(!isnan(Δs))
    @assert(!isnan(initial_speed))

    Δt = Δs / initial_speed

    if isnan(Δt)
        println("initial_carind: ", initial_carind)
        println("velFx: ", velFx)
        println("velFy: ", velFy)
        println("seg: ", seg)
        println("Δs: ", Δs)
        println("Δd: ", Δd)
        println("initial_speed: ", initial_speed)
        println(project_point_to_streetmap(posGx_A, posGy_A, sn))
        println(project_point_to_streetmap(posGx_B, posGy_B, sn))
    end
    @assert(!isnan(Δt))

    (Δt, Δd)
end

function allocate_empty_histobin(histobin_params::ParamsHistobin)
    n_bins_x = nlabels(histobin_params.discx)
    n_bins_y = nlabels(histobin_params.discy)
    zeros(Float64, n_bins_x, n_bins_y)
end
function update_histobin!(
    histobin::Matrix{Float64},
    pdset::PrimaryDataset,
    sn::StreetNetwork,
    seg::PdsetSegment,
    histobin_params::ParamsHistobin,
    )

    Δt, Δy = calc_trace_deviation(pdset, sn, seg)

    bin_x = encode(histobin_params.discx, Δt)
    bin_y = encode(histobin_params.discy, Δy)

    histobin[bin_x, bin_y] += 1.0

    histobin
end
function calc_histobin(
    pdsets::Vector{PrimaryDataset},
    streetnets::Vector{StreetNetwork},
    pdset_segments::Vector{PdsetSegment},
    histobin_params::ParamsHistobin,
    foldset::FoldSet, # :seg
    )

    histobin = allocate_empty_histobin(histobin_params)

    for index in foldset
        seg = pdset_segments[index]
        update_histobin!(histobin,
                         pdsets[seg.pdset_id],
                         streetnets[seg.streetnet_id],
                         seg, histobin_params)
    end

    histobin
end
function calc_histobin(
    pdsets::Vector{PrimaryDataset},
    streetnets::Vector{StreetNetwork},
    pdset_segments::Vector{PdsetSegment},
    histobin_params::ParamsHistobin
    )

    histobin = allocate_empty_histobin(histobin_params)

    for seg in pdset_segments
        update_histobin!(histobin,
                         pdsets[seg.pdset_id],
                         streetnets[seg.streetnet_id],
                         seg, histobin_params)
    end

    histobin
end

function calc_histobin(Δt_arr::Vector{Float64}, Δy_arr::Vector{Float64}, histobin_params::ParamsHistobin)

    histobin = allocate_empty_histobin(histobin_params)

    for (Δt, Δy) in zip(Δt_arr, Δy_arr)

        bin_x = encode(histobin_params.discx, Δt)
        bin_y = encode(histobin_params.discy, Δy)

        histobin[bin_x, bin_y] += 1.0
    end

    histobin
end

function KL_divergence_categorical(histobinA::Matrix{Float64}, histobinB::Matrix{Float64})

    Ap = histobinA ./ sum(histobinA)
    Bp = histobinB ./ sum(histobinB)

    KL_divergence = 0.0
    for i = 1 : length(Ap)
        KL_divergence += Ap[i]*log(Ap[i]/Bp[i])
    end
    KL_divergence::Float64
end
function KL_divergence_dirichlet(histobinA::Matrix{Float64}, histobinB::Matrix{Float64})
    α0 = sum(histobinA)
    KL_divergence = 0.0
    KL_divergence += lgamma(α0)
    KL_divergence -= lgamma(sum(histobinB))
    KL_divergence -= sum([lgamma(a) for a in histobinA])
    KL_divergence += sum([lgamma(b) for b in histobinB])
    for i = 1 : length(histobinA)
        KL_divergence += (histobinA[i] - histobinB[i])*(digamma(histobinA[i]) - digamma(α0))
    end
    KL_divergence::Float64
end
function KL_divergence_univariate_gaussian(μ1::Float64, μ2::Float64, σ1::Float64, σ2::Float64)
    Δμ = (μ1-μ2)
    log(σ2/σ1) + (σ1*σ1 + Δμ*Δμ)/(2σ2*σ2) - 0.5
end