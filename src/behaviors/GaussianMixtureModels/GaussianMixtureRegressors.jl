module GaussianMixtureRegressors

# http://cs229.stanford.edu/section/more_on_gaussians.pdf
# http://pypr.sourceforge.net/mog.html

using ProbDrive2016ITSC
using GaussianMixtures
using Distributions
using PyCall

@pyimport sklearn.mixture as mixture

import ProbDrive2016ITSC:
    AbstractVehicleBehavior,
    AbstractVehicleBehaviorPreallocatedData,
    AbstractVehicleBehaviorTrainParams,

    preallocate_learning_data,
    select_action,
    calc_action_loglikelihood,
    train

export
    GMRBehavior,
    GMR_TrainParams,

    select_action,
    calc_action_loglikelihood,
    train

const DEFAULT_MIN_EIGENVALUE = 1e-6

function ensure_posdef!(Σ::Matrix{Float64}, min_eigenvalue::Float64)

    @assert(size(Σ, 1) == size(Σ, 2))

    Λ, V = eig(Σ)
    for (j,λ) in enumerate(Λ)
        Λ[j] = max(DEFAULT_MIN_EIGENVALUE, λ)
    end
    Σ[:] = V*diagm(Λ)/V

    Σ
end

function nearestSPD(A::Matrix{Float64})

    # see http://www.mathworks.com/matlabcentral/fileexchange/42885-nearestspd

    # output:
    #  α, β ≥ 0.0 such that
    #         α ≤ δ₂(A) ≤ β ≤ α + 2 max(fα, tol)
    #  and a PSD matrix X such that |A - X|₂ = β

    n = size(A, 1)
    @assert(n == size(A, 2)) # ensure it is square

    I = eye(n)

    # symmetrize A into B
    B = (A+A')./2

    # Compute the symmetric polar factor of B. Call it H.
    # Clearly H is itself SPD.
    U, σ, V = svd(B)
    H = V*diagm(σ)*V'

    # get Ahat in the above formula
    Ahat = (B+H)/2

    # ensure symmetry
    Ahat = (Ahat + Ahat')/2;

    # test that Ahat is in fact PD. if it is not so, then tweak it just a bit.
    worked = false
    iteration_count = 0
    while !worked && iteration_count < 100

        iteration_count += 1

        try
            chol(Ahat)
            worked = true
        catch
            # do nothing
        end

        if !worked
            # Ahat failed the chol test. It must have been just a hair off,
            # due to floating point trash, so it is simplest now just to
            # tweak by adding a tiny multiple of an identity matrix.

            min_eig = minimum(eigvals(Ahat))
            Ahat = Ahat + (-min_eig*iteration_count.^2 + eps(Float32))*I
        end
    end

    Ahat
end

type GMR_TrainParams <: AbstractVehicleBehaviorTrainParams

    targets::ModelTargets
    indicators::Vector{AbstractFeature}

    # params for sci-kit learn GMM fitting
    n_components::Int
    n_gmm_iter::Int # number of iterations taken during GMM fitting
    n_init::Int # number of GMMs to attempt to fit
    verbose::Int # defaults to 0, with 1 prints some stuff, with 2 prints more
    tol::Float64 # tolerance on convergence, if Δlogl < tol, we have converged
    min_covar::Float64 # minimum size of covariance matrix diagonal entry

    # params for julia fitting and training
    max_n_indicators::Int
    use_PCA::Bool
    unlearned_component_weight::Float64 # the weighting given to the unlearned component

    function GMR_TrainParams(;
        targets::ModelTargets = ModelTargets(FUTUREDESIREDANGLE, FUTUREACCELERATION),
        indicators::Vector{AbstractFeature} = AbstractFeature[],

        # params for sci-kit learn GMM fitting
        n_components::Int = 2,
        n_gmm_iter::Int = 100,
        n_init::Int = 3,
        verbose::Int = 0,
        tol::Float64 = 1e-3,
        min_covar::Float64 = 1e-6,

        # params for julia fitting and training
        max_n_indicators::Int = 3,
        use_PCA::Bool = false,
        unlearned_component_weight::Float64 = 0.01, # ∈ [0,1]
        )

        retval = new()

        retval.targets = targets
        retval.indicators = indicators

        retval.n_components = n_components
        retval.n_gmm_iter = n_gmm_iter
        retval.n_init = n_init
        retval.verbose = verbose
        retval.tol = tol
        retval.min_covar = min_covar

        retval.max_n_indicators = max_n_indicators
        retval.use_PCA = use_PCA
        retval.unlearned_component_weight = unlearned_component_weight

        retval
    end
end
function Base.print(io::IO, θ::GMR_TrainParams)
    println(io, "GMR Train Params")
    println(io, "\ttargets: ", θ.targets)
    println(io, "\tindicators: ", map(f->symbol(f), θ.indicators))
    println(io, "\tn_components:     ", θ.n_components)
    println(io, "\tmax_n_indicators: ", θ.max_n_indicators)
    println(io, "\tuse_PCA:   ", θ.use_PCA)
    println(io, "\tmin_covar:        ", θ.min_covar)
end

type GMR_PreallocatedData <: AbstractVehicleBehaviorPreallocatedData

    YX::Matrix{Float64} # [nframes × 2+p]
                        # this includes both targets as 1st two columns (lat, lon)
                        # and then the predictors
    YZ::Matrix{Float64} # [nframes × 2+p]
                        # this uses PCA'd predictors

    extractor::FeatureSubsetExtractor
    preprocess::DataPreprocessor
    preprocess_pca::DataPreprocessor
    preprocess_target::DataPreprocessor

    function GMR_PreallocatedData(dset::ModelTrainingData2, params::GMR_TrainParams)

        retval = new()

        targets = params.targets
        indicators = params.indicators
        trainingframes = dset.dataframe
        nframes = nrow(trainingframes)
        nindicators = length(indicators)


        X = Array(Float64, nindicators, nframes)
        Y = Array(Float64, 2, nframes)

        pull_design_and_target_matrices!(X, Y, trainingframes, targets, indicators)

        ###########################

        preprocess_target = ChainedDataProcessor(2)
        push!(preprocess_target, Y, DataClamper)
        process!(Y, preprocess_target.processors[end]) # process in place
        # push!(preprocess_target, Y, DataScaler) # ensure data is standardized before running PCA
        # process!(Y, preprocess_target.processors[end]) # process in place

        @assert(findfirst(v->isnan(v), Y) == 0)
        @assert(findfirst(v->isinf(v), Y) == 0)

        extractor = FeatureSubsetExtractor(Array(Float64, nindicators), indicators)
        preprocess = ChainedDataProcessor(extractor)
        push!(preprocess, X, DataNaReplacer, extractor.indicators)
        process!(X, preprocess.processors[end]) # process in place

        @assert(findfirst(v->isnan(v), X) == 0)
        @assert(findfirst(v->isinf(v), X) == 0)

        push!(preprocess, X, DataClamper)
        process!(X, preprocess.processors[end]) # process in place

        @assert(findfirst(v->isnan(v), X) == 0)
        @assert(findfirst(v->isinf(v), X) == 0)

        push!(preprocess, X, DataScaler) # ensure data is standardized before running PCA
        process!(X, preprocess.processors[end]) # process in place

        @assert(findfirst(v->isnan(v), X) == 0)
        @assert(findfirst(v->isinf(v), X) == 0)

        retval.YX = vcat(Y, copy(X))
        @assert(findfirst(v->isnan(v), retval.YX) == 0)
        @assert(findfirst(v->isinf(v), retval.YX) == 0)

        # Add a PCA scaler
        preprocess_pca = deepcopy(preprocess)
        n_PCA_features = min(nindicators, 20)
        push!(preprocess_pca, X, DataLinearTransform, n_PCA_features)

        Z = Array(Float64, n_PCA_features, nframes)
        Z = process!(Z, X, preprocess_pca.processors[end])

        retval.YZ = vcat(Y, Z)
        @assert(findfirst(v->isnan(v), retval.YZ) == 0)
        @assert(findfirst(v->isinf(v), retval.YZ) == 0)

        retval.extractor = extractor
        retval.preprocess = preprocess
        retval.preprocess_pca = preprocess_pca
        retval.preprocess_target = preprocess_target

        retval
    end
end
function preallocate_learning_data(
    dset::ModelTrainingData2,
    params::GMR_TrainParams)

    GMR_PreallocatedData(dset, params)
end

######################################################################
#
#
#
######################################################################

type GMR
    converged::Bool

    # μ₁₋₂ = μ₁ + Σ₁₂ * Σ₂₂⁻¹ * (x₂ - μ₂) = A*x₂ + b
    vec_A::Vector{Matrix{Float64}} # [n_components [2×nindicators]]
    vec_b::Vector{Vector{Float64}} # [n_components [2]]

    # pdf(a|p), means μⱼ_ₚ and weights βⱼ(p) are functions of p and must be updated every time, covariance is constant
    # the last component is the extra, unlearned component, and its weight is static
    mixture_Act_given_Obs::MixtureModel{Multivariate,Continuous,MvNormal}

    # pdf(p), all pre-computed. Used to compute βⱼ(p)
    mixture_Obs::MixtureModel{Multivariate,Continuous,MvNormal} # p(obs), all pre-computed, should never be edited

    function GMR(gmm::PyObject, params::GMR_TrainParams)

        retval = new()

        if !gmm[:converged_]
            retval.converged = false
        else

            weights = gmm[:weights_]::Vector{Float64} # [n_components]
            means = gmm[:means_]::Matrix{Float64}     # [n_components, n_features]
            covars = gmm[:covars_]::Array{Float64, 3} # shape depends on covariance_type
            @assert(gmm[:covariance_type] == "full")
            # (n_components, n_features)             if 'spherical',
            # (n_features, n_features)               if 'tied',
            # (n_components, n_features)             if 'diag',
            # (n_components, n_features, n_features) if 'full'

            @assert(0.0 ≤ params.unlearned_component_weight ≤ 1.0)

            n_targets = 2
            n_indicators = size(means, 2) - n_targets
            n_components = length(weights)

            # println("n_targets:    ", n_targets)
            # println("n_indicators: ", n_indicators)
            # println("n_components: ", n_components)

            vec_A = Array(Matrix{Float64}, n_components) # μ₁₋₂ = μ₁ + Σ₁₂ * Σ₂₂⁻¹ * (x₂ - μ₂) = A*x₂ + b
            vec_b = Array(Vector{Float64}, n_components)
            vec_G = Array(MvNormal, n_components)
            vec_H = Array(MvNormal, n_components)

            for i = 1 : n_components

                μₐ = vec(means[i,1:n_targets])
                μₚ = vec(means[i,n_targets+1:end])

                Σ = reshape(covars[i, :, :], (n_targets + n_indicators,n_targets + n_indicators))

                Σₐₐ = Σ[1:n_targets,1:n_targets]
                Σₐₚ = Σ[1:n_targets,n_targets+1:end]
                Σₚₚ = nearestSPD(Σ[n_targets+1:end,n_targets+1:end])
                iΣₚₚ = inv(Σₚₚ)

                A = Σₐₚ * iΣₚₚ
                vec_A[i] = A
                vec_b[i] = vec(μₐ - A*μₚ)
                C = nearestSPD(Σₐₐ - Σₐₚ * iΣₚₚ * ((Σₐₚ)'))

                vec_G[i] = MvNormal(Array(Float64, n_targets), C) # p(action|obs), mean and weighting must be updated with each observation, cov is pre-computed
                vec_H[i] = MvNormal(μₚ, Σₚₚ) # p(obs), all pre-computed, should never be edited
            end

            # add the un-learned component (needed because direct learning will NaN with outliers)
            push!(vec_G, MvNormal([0.0,0.0], [0.01,0.1])) # TODO(tim): tune this
            wᵤ = params.unlearned_component_weight
            Act_given_Obs_prior = push!(fill((1-wᵤ)/n_components, n_components), wᵤ)

            retval.converged = true
            retval.vec_A = vec_A
            retval.vec_b = vec_b
            retval.mixture_Act_given_Obs = MixtureModel(vec_G, Act_given_Obs_prior) # p(action|obs), mean and weighting must be updated with each observation, cov is pre-computed
            retval.mixture_Obs = MixtureModel(vec_H, weights) # p(obs), all pre-computed, should never be edited
        end

        retval
    end
end
function Base.print(model::GMR)
    if !model.converged
        println("GMR: <not converged>")
    else
        println("GMR:")
        for (i, mat) in enumerate(model.vec_A)
            println(i)
            print("\t[")
            for j in 1:size(mat,2)
                @printf("  %10.6f", mat[1,j])
            end
            @printf("] + [  %10.6f]\n", model.vec_b[i][1])
            print("\t[")
            for j in 1:size(mat,2)
                @printf("  %10.6f", mat[2,j])
            end
            @printf("] + [  %10.6f]\n", model.vec_b[i][2])
        end
        println("\tmixture_Obs: ")
        println("\t\tprior: ", model.mixture_Obs.prior)
    end
end

_n_learned_components(gmr::GMR) = length(gmr.vec_A)
function nsuffstats(gmr::GMR)


    dimA = length(gmr.vec_A[1])
    n_learned_components = _n_learned_components(gmr)

    n_learned_components * (2*dimA + 2 # bias
                                   + 3 # covariance in mixture_Act_given_Obs
                                   + div(dimA*dimA,2)) # covariance for mixture_Obs
end
function _calc_action_loglikelihood(gmr::GMR, a::Vector{Float64}, f::Vector{Float64})

    if !gmr.converged
        println("converged: ", gmr.converged)
    end

    mixture_Act_given_Obs = gmr.mixture_Act_given_Obs
    mixture_Obs = gmr.mixture_Obs
    n_learned_components = _n_learned_components(gmr)

    for j in 1 : n_learned_components

        # compute the β value, unweighted
        # βⱼ(f) ∝ wⱼ Nⱼ(μₚ, Σₚ)
        wⱼ = mixture_Obs.prior.p[j]
        Nⱼ = mixture_Obs.components[j]

        # println("a   ", a)
        # println("f   ", f)
        # println("wⱼ   ", wⱼ)
        # println("Nⱼ\n", Nⱼ)

        mixture_Act_given_Obs.prior.p[j] = wⱼ * pdf(Nⱼ, f)

        # compute the conditional mean
        #  μₐ_ₚ = A⋅f + b
        A = gmr.vec_A[j]
        b = gmr.vec_b[j]

        # println("A:\n", A)
        # println("b:\n", b)

        mixture_Act_given_Obs.components[j].μ[:] = A*f + b
    end

    # normalize the β values
    sum_β = sum(mixture_Act_given_Obs.prior.p)
    for i in 1 : n_learned_components
        mixture_Act_given_Obs.prior.p[i] /= sum_β
    end

    if findfirst(v->isinf(v), mixture_Act_given_Obs.prior.p) != 0 || findfirst(v->isnan(v), mixture_Act_given_Obs.prior.p) != 0

        @printf("j  %10s  %10s\n", "wⱼ", "pdf(Nⱼ, f)")
        for j in 1 : n_learned_components
            wⱼ = mixture_Obs.prior.p[j]
            Nⱼ = mixture_Obs.components[j]
            @printf("%2d  %10.6f  %10.6f\n", j, wⱼ, pdf(Nⱼ, f))
        end
        # sleep(100)

        println("INF1: ")
        println("\t a|p prior:  ", mixture_Act_given_Obs.prior.p)
        println("\t sum_β:      ", sum_β)
        println("\t indicators: ", f)
        println("\t targets:    ",    a)

        sleep(0.01)

        return -Inf
    end

    # compute logl using the built-in mixture
    logl = logpdf(mixture_Act_given_Obs, a)

    if isinf(logl) || isnan(logl)
        println(mixture_Act_given_Obs.prior.p)
        println(sum_β)
        println("INF2: obs: ", f, "  action: ",  a)
        sleep(0.01)
    end

    logl
end
function _train_model(gmm::PyObject, YX::Matrix{Float64}, chosen_indicators::Vector{Int}, params::GMR_TrainParams)
    # println("GMM trainging ", chosen_indicators)
    # tic()
    gmm[:fit](YX[:,[1;2;chosen_indicators+2]])
    # toc()
    GMR(gmm, params)
end
function _calc_bic_score(gmr::GMR, YX::Matrix{Float64}, chosen_indicators::Vector{Int})

    bic = -Inf

    if gmr.converged

        a = Array(Float64, 2)
        f = Array(Float64, length(chosen_indicators))

        m = size(YX, 1)
        logl = 0.0
        for i in 1 : m

            a[1] = YX[i,1]
            a[2] = YX[i,2]

            for (j,p) in enumerate(chosen_indicators)
                f[j] = YX[i, 2+p]
            end

            logl += _calc_action_loglikelihood(gmr, a, f)
        end

        # # DEBUG
        # println("\tlogl:     ", logl)
        # println("\t2nd term: ", log(m)*nsuffstats(gmr)/2)
        # println("\ttotal:    ", logl - log(m)*nsuffstats(gmr)/2)

        bic = logl - log(m)*nsuffstats(gmr)/2
    end

    bic
end

######################################################################
#
#
#
######################################################################

type GMRBehavior <: AbstractVehicleBehavior

    #=
    The GMR behavior does Gaussian Mixture Regression

    1 - fit a joint GMM over (lat,lon,p₁,p₁,...,pₙ)

        pdf({a,P}) = pdf(x) = ∑ wⱼ N(μⱼ,Σⱼ)

    2 - convert to a conditional: (lat,lon | p₁,p₁,...,pₙ)

        pdf(a|P) = ∑ βⱼ(p) N(μⱼ_ₚ, Σⱼ_ₚ)

            μⱼ_ₚ = μₐ + Σₐₚ/Σₚₚ (p - μₚ)
            Σⱼ_ₚ = Σₐₐ - (Σₐₚ/Σₚₚ)Σₚₐ
            βⱼ(p) ∝ wⱼ N(μₚ, Σₚₚ)
    =#

    gmr::GMR
    targets::ModelTargets
    extractor::FeatureSubsetExtractor # - note that the extractor and processor are linked internally
    processor::DataPreprocessor
    processor_target::DataPreprocessor

    function GMRBehavior(
        gmr::GMR,
        chosen_indicators::Vector{Int},
        params::GMR_TrainParams,
        pre::GMR_PreallocatedData,
        )

        indicators = pre.extractor.indicators[chosen_indicators]
        preprocess = params.use_PCA > 0 ? pre.preprocess_pca : pre.preprocess

        retval = new()
        retval.gmr = gmr
        retval.targets = params.targets
        retval.extractor = FeatureSubsetExtractor(deepcopy(pre.extractor.indicators))
        retval.processor = deepcopy(preprocess, retval.extractor)
        push!(retval.processor, DataSubset, chosen_indicators)
        retval.processor_target = deepcopy(pre.preprocess_target)
        retval
    end
end
function Base.print(io::IO, GM::GMRBehavior)

    num_components = ncomponents(GM.gmr.mixture_Obs)
    n_indicators = length(GM.processor.z)

    println(io, "Gaussian Mixture Regression Model")
    println(io, "\ttargets: ", GM.targets)
    println(io, "\tncomponents: ", num_components)
    println(io, "\tnindicators: ", n_indicators)
    print(GM.gmr)
end
ProbDrive2016ITSC.get_targets(behavior::GMRBehavior) = behavior.targets
function ProbDrive2016ITSC.get_indicators(behavior::GMRBehavior)
    subset_index = findfirst(p->isa(p, DataSubset), behavior.processor.processors)
    subset = behavior.processor.processors[subset_index]
    behavior.extractor.indicators[subset.indeces]
end
_n_learned_components(GM::GMRBehavior) = length(GM.gmr.vec_A)

######################################################################
#
#
#
######################################################################

const TOLERANCE_ZERO_PROB_WEIGHT = 1e-20

function select_action(
    behavior::GMRBehavior,
    runlog::RunLog,
    sn::StreetNetwork,
    colset::UInt,
    frame::Int
    )

    observe!(behavior.extractor, runlog, sn, colset, frame)
    process!(behavior.processor)

    f = behavior.processor.z

    # condition on the mixture
    vec_A = behavior.gmr.vec_A
    vec_b = behavior.gmr.vec_b
    n_components = _n_learned_components(behavior)
    mixture_Obs = behavior.gmr.mixture_Obs
    mixture_Act_given_Obs = behavior.gmr.mixture_Act_given_Obs

    total_prior_sum = mixture_Act_given_Obs.prior.p[end] # include weight from unlearned component
    for i = 1 : n_components
        mixture_Act_given_Obs.prior.p[i] = mixture_Obs.prior.p[i] * pdf(mixture_Obs.components[i], f)
        total_prior_sum += mixture_Act_given_Obs.prior.p[i]
    end
    if abs(total_prior_sum) < TOLERANCE_ZERO_PROB_WEIGHT
        return (0.0,0.0)
    end

    for i = 1 : n_components
        mixture_Act_given_Obs.prior.p[i] /= total_prior_sum
    end
    @assert(!isnan(mixture_Act_given_Obs.prior.p[1]))


    # pick a component to sample from
    component_index = rand(mixture_Act_given_Obs.prior)
    component = mixture_Act_given_Obs.components[component_index]

    if component_index != n_components+1
        # not the unlearned component
        # modify said component with observation
        component.μ[:] = vec_b[component_index] + vec_A[component_index]*f
    end

    # sample from said component
    Distributions._rand!(component, behavior.processor_target.x)
    process!(behavior.processor_target)

    # return the selected action
    action_lat = behavior.processor_target.z[1]
    action_lon = behavior.processor_target.z[2]

    (action_lat, action_lon)
end

function _set_and_process_action!(behavior::GMRBehavior, action_lat::Float64, action_lon::Float64)
    behavior.processor_target.x[1] = action_lat
    behavior.processor_target.x[2] = action_lon
    process!(behavior.processor_target)
    behavior.processor_target
end
function _calc_action_loglikelihood(behavior::GMRBehavior)

    # NOTE(tim): this assumes that:
    #   - observations have already been pulled and are in behavior.processor.z
    #   - action_lat and action_lon are already processed and in behavior.processor_target.z

    a = behavior.processor_target.z
    f = behavior.processor.z
    _calc_action_loglikelihood(behavior.gmr, a, f)
end

function calc_action_loglikelihood(
    behavior::GMRBehavior,
    runlog::RunLog,
    sn::StreetNetwork,
    colset::UInt,
    frame::Int,
    action_lat::Float64,
    action_lon::Float64,
    )

    #=
    Compute the log-likelihood of the action taken during a single frame
    given the VehicleBehaviorGaussian.
    =#

    observe!(behavior.extractor, runlog, sn, colset, frame)
    process!(behavior.processor)

    _set_and_process_action!(behavior, action_lat, action_lon)

    _calc_action_loglikelihood(behavior)
end
function calc_action_loglikelihood(
    behavior::GMRBehavior,
    features::DataFrame,
    frameind::Integer,
    )

    observe!(behavior.extractor, features, frameind)
    process!(behavior.processor)

    action_lat = features[frameind, symbol(behavior.targets.lat)]::Float64
    action_lon = features[frameind, symbol(behavior.targets.lon)]::Float64
    # println("action_lat: ", action_lat)
    # println("action_lon: ", action_lon)
    _set_and_process_action!(behavior, action_lat, action_lon)

    # println("after: behavior.processor_target.x: ", behavior.processor_target.x)
    # println("after: behavior.processor_target.z: ", behavior.processor_target.z)

    _calc_action_loglikelihood(behavior)
end

######################################################################
#
#
#
######################################################################

function train(
    training_data::ModelTrainingData2, # TODO: do we even need this?
    preallocated_data::GMR_PreallocatedData,
    params::GMR_TrainParams,
    foldset::FoldSet,
    )

    YX = params.use_PCA ?
            copy_matrix_fold(preallocated_data.YZ, foldset) :
            copy_matrix_fold(preallocated_data.YX, foldset)

    @assert(findfirst(v->isnan(v), YX) == 0)
    @assert(findfirst(v->isinf(v), YX) == 0)

    # -----------------------------------------
    # run greedy ascent
    #  - always add the next best parent

    gmm = mixture.GMM(n_components    = params.n_components, # int, optional. Number of mixture components. Defaults to 1.
                      covariance_type = "full", # string, optional. ‘spherical’, ‘tied’, ‘diag’, ‘full’. Defaults to ‘diag’
                      random_state    = 1, # RandomState or an int seed (None by default)
                      min_covar       = params.min_covar, # float, optional
                      n_iter          = params.n_gmm_iter,
                      n_init          = params.n_init,
                      tol             = params.tol,
                      verbose         = params.verbose
                      )

    n_indicators = size(YX, 2) - 2
    max_n_indicators = min(params.max_n_indicators, n_indicators, div(size(YX, 1), 2))
    chosen_indicators = Int[] # start with no parents
    best_model = _train_model(gmm, YX, chosen_indicators, params)
    best_score = _calc_bic_score(best_model, YX, chosen_indicators)

    finished = false
    while !finished && length(chosen_indicators) < max_n_indicators

        finished = true
        test_indicators = [deepcopy(chosen_indicators); 0]

        for i in 1:n_indicators
            if !in(i, chosen_indicators)
                test_indicators[end] = i
                test_model = _train_model(gmm, YX, test_indicators, params)
                test_score = _calc_bic_score(test_model, YX, test_indicators)

                if test_score > best_score
                    best_score = test_score
                    best_model = test_model
                    finished = false
                    chosen_indicators = deepcopy(test_indicators)
                end
            end
        end
    end

    GMRBehavior(best_model, chosen_indicators, params, preallocated_data)
end

end # end module