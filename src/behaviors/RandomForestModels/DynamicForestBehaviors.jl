module DynamicForestBehaviors

using ProbDrive2016ITSC

using RandomForestBehaviors.MvDecisionTrees
using Distributions.MvNormal

import ProbDrive2016ITSC:
    AbstractVehicleBehavior,
    AbstractVehicleBehaviorTrainParams,
    AbstractVehicleBehaviorPreallocatedData,

    preallocate_learning_data,
    select_action,
    calc_action_loglikelihood,
    train

export
    DynamicForestBehavior,
    DF_TrainParams,
    AutoregressiveMvNormLeaf,

    select_action,
    calc_action_loglikelihood,
    train

global GLOBAL_DIAG_SHRINKAGE_THRESOLD = 1e-6
DEFAULT_AUTOREGRESSION_CONSTANT = 0.5
DEFAULT_NUM_AUTOREGRESSION_PREDICTORS = 3
DEFAULT_NUM_PREDICTOR_SAMPLES = 10

################################################################

immutable AutoregressiveMvNormLeaf <: MvDecisionTrees.Leaf
    A::Matrix{Float64} # regression matrix A = UΦᵀ(ΦΦᵀ + γI)⁻¹ where γ is a ridge regularization param, U is column-wise targets, Φ is column-wise features
    m::MvNormal        # μ = Aϕ must be overwritten on each use
    predictor_indeces::Vector{Int} # indeces of predictors
end

function _regress_on_predictor_subset{T<:AbstractFloat, S<:Real}(
    data::MvDecisionTrees.TreeData{T,S},
    assignment_id::Int,
    Γ::Matrix{T}, # ϕ+1 × ϕ+1 diagonal γ matrix (preallocated)
    U::Matrix{T}, #   o × m
    Φ::Matrix{T}, # ϕ+1 × m (is overwritten)
    A::Matrix{T}, #   o × ϕ+1 (is overwritten)
    predictor_indeces::Vector{Int}, # length ϕ
    )

    #=
    Returns the loss when using this predictor set
    Overwrites Φ
    =#

    n = size(data.X,1)
    o, m = size(U)
    ϕ = length(predictor_indeces)
    @assert(size(Φ,1) == ϕ+1)
    @assert(size(Φ,2) == m)
    @assert(size(Γ,1) == size(Γ,2) == ϕ+1 )

    u_index = 0
    for i in 1 : n
        if data.assignment[i] == assignment_id
            u_index += 1
            for (j,k) in enumerate(predictor_indeces)
                Φ[j,u_index] = data.X[i,k] # features : [nrows, nfeatures], Φ: [nfeatures, nrows]
            end
            Φ[ϕ+1,u_index] = 1.0 # constant 1-predictor
        end
    end

    # regress

    A[:] = U*Φ'/(Φ*Φ' + Γ)

    # calc residual loss (square deviation)

    loss = 0.0
    for i in 1 : n
        if data.assignment[i] == assignment_id
            for j in 1 : o
                r = data.Y[j,i] - A[j,ϕ+1]
                for k in 1 : ϕ
                    r -= A[j,k]*data.X[i,predictor_indeces[k]]
                end
                loss += r*r # squared error
            end
        end
    end

    loss
end
function MvDecisionTrees.build_leaf{T<:AbstractFloat, S<:Real}(::Type{AutoregressiveMvNormLeaf},
    data::MvDecisionTrees.TreeData{T,S},
    assignment_id::Int,
    γ::T=DEFAULT_AUTOREGRESSION_CONSTANT, # ridge-regularization parameter
    n_autoregression_predictors::Int = DEFAULT_NUM_AUTOREGRESSION_PREDICTORS,
    n_random_predictor_samples::Int = DEFAULT_NUM_PREDICTOR_SAMPLES,
    )

    #=
    Fit the regression matrix A = UΦᵀ(ΦΦᵀ + γI)⁻¹
       γ is a ridge regularization param ≥ 0
         0 → least squares
         large γ → increased penalty on regression coefficients
       U is column-wise targets
       Φ is column-wise features plus a constant 1-feature

    This function picks a set of 'n_autoregression_predictors' predictors used to regress on the targets
    This is done 'n_random_predictor_samples' times
    The best in terms of residual error is kept
    The covariance matrix is static and computed based on the final residual errors
    =#

    @assert(n_autoregression_predictors ≥ 0)
    @assert(n_random_predictor_samples > 0)

    n, o, p = size(data)
    id_count = MvDecisionTrees.get_id_count(data, assignment_id)
    @assert(id_count > 0)

    ##################################
    # pull target matrix

    U = Array(T, o, id_count) # column-wise contatenation of target variables

    u_index = 0
    for i in 1 : n
        if data.assignment[i] == assignment_id
            u_index += 1
            for j in 1 : o
                U[j,u_index] = data.Y[j,i]
            end
        end
    end

    #################################
    # pull predictor matrix

    ϕ = min(n_autoregression_predictors, p) # number of real features
    Φ = Array(T, ϕ+1, id_count) # column-wise contatenation of predictors (also including constant 1-predictor)
    Γ = diagm(fill(γ, (ϕ+1)))
    A = Array(T, o, ϕ+1)
    predictor_indeces = Array(Int, ϕ)

    if ϕ == 0
        _regress_on_predictor_subset(data, assignment_id, Γ, U, Φ, A, predictor_indeces) # just solve for A
    else
        # pick several samples of predictor_indeces
        temp_A = deepcopy(A)
        temp_predictor_indeces = deepcopy(predictor_indeces)
        best_loss = Inf

        for i = 1 : n_random_predictor_samples
            MvDecisionTrees._reservoir_sample!(temp_predictor_indeces, p, ϕ)
            loss = _regress_on_predictor_subset(data, assignment_id, Γ, U, Φ, temp_A, temp_predictor_indeces)
            if loss < best_loss
                best_loss = loss
                copy!(A, temp_A)
                copy!(predictor_indeces, temp_predictor_indeces)
            end
        end

        # in case we found nothing
        if isinf(best_loss)
            fill!(A, 0.0)
            fill!(predictor_indeces, 1)
        end
    end

    ###########################################################3
    # calc covariance

    Σ = zeros(T, o, o)
    r_vec = Array(Float64, o)

    # obtain upper triangle
    for i in 1 : n
        if data.assignment[i] == assignment_id

            for j in 1 : o
                r_vec[j] = data.Y[j,i] - A[j,ϕ+1]
                for k in 1 : ϕ
                    r_vec[j] -= A[j,k]*data.X[i,predictor_indeces[k]]
                end
            end

            for a in 1 : o
                for b in a : o
                    Σ[a,b] += r_vec[a] * r_vec[b]
                end
            end
        end
    end

    if id_count < 2
        warn("Not enough samples to build a covariance matrix")
        for i in 1 : o
            Σ[i,i] = 1.0
        end
    else
        for i in 1 : o*o
            @assert(!isnan(Σ[i]))
            Σ[i] /= (id_count-1)
        end

        ε = GLOBAL_DIAG_SHRINKAGE_THRESOLD
        MvDecisionTrees._diagonal_shrinkage!(Σ, ε)
        for a in 2:o
            for b in 1:a-1
                Σ[a,b] = Σ[b,a]
            end
        end
    end

    AutoregressiveMvNormLeaf(A, MvNormal(Array(T, o), Σ), predictor_indeces)
end

################################################################

type DF_TrainParams <: AbstractVehicleBehaviorTrainParams

    targets::ModelTargets
    indicators::Vector{AbstractFeature}

    ntrees::Int
    max_tree_depth::Int
    n_split_tries::Int
    min_samples_split::Int
    min_samples_leaves::Int
    n_random_predictor_samples::Int
    n_autoregression_predictors::Int
    use_PCA::Bool

    diag_shrinkage_threshold::Float64
    min_split_improvement::Float64
    partial_sampling::Float64
    autogression_coef::Float64

    function DF_TrainParams(;
        targets::ModelTargets = ModelTargets(Features.FUTUREDESIREDANGLE, Features.FUTUREACCELERATION),
        indicators::Vector{AbstractFeature} = AbstractFeature[],

        ntrees::Integer=10,
        max_tree_depth::Integer=5,
        n_split_tries::Integer=50,
        min_samples_split::Integer=100,
        min_samples_leaves::Integer=20,
        n_random_predictor_samples::Integer=10,
        n_autoregression_predictors::Integer=2,
        use_PCA::Bool=false,

        diag_shrinkage_threshold::Float64=1e-6,
        min_split_improvement::Float64=0.0,
        partial_sampling::Float64=0.7,
        autogression_coef::Float64=0.1,
        )

        retval = new()
        retval.targets = targets
        retval.indicators = indicators
        retval.ntrees = ntrees
        retval.max_tree_depth = max_tree_depth
        retval.n_split_tries = n_split_tries
        retval.min_samples_split = min_samples_split
        retval.min_samples_leaves = min_samples_leaves
        retval.n_random_predictor_samples = n_random_predictor_samples
        retval.n_autoregression_predictors = n_autoregression_predictors
        retval.use_PCA = use_PCA

        retval.diag_shrinkage_threshold = diag_shrinkage_threshold
        retval.min_split_improvement = min_split_improvement
        retval.partial_sampling = partial_sampling
        retval.autogression_coef = autogression_coef
        retval
    end
end
type DF_PreallocatedData <: AbstractVehicleBehaviorPreallocatedData

    # the full dataset matrices pre-extracted
    X::Matrix{Float64} # [nfeatures × nframes]
    Y::Matrix{Float64} # [2 × nframes]
    Z::Matrix{Float64} # [n_pca_features × nframes]

    extractor::FeatureSubsetExtractor
    preprocess::DataPreprocessor
    preprocess_pca::DataPreprocessor
    action_clamper::DataClamper

    function DF_PreallocatedData(dset::ModelTrainingData2, params::DF_TrainParams)

        retval = new()

        targets = params.targets
        indicators = params.indicators
        trainingframes = dset.dataframe
        nframes = nrow(trainingframes)
        nindicators = length(indicators)

        X = Array(Float64, nindicators, nframes) # NOTE(tim): initially indicators × nframes, we later transpose it
        Y = Array(Float64, 2, nframes)

        pull_design_and_target_matrices!(X, Y, trainingframes, targets, indicators)

        ######################

        extractor = FeatureSubsetExtractor(indicators)
        preprocess = ChainedDataProcessor(extractor)
        push!(preprocess, X, DataNaReplacer, extractor.indicators)
        process!(X, preprocess.processors[end]) # process in place
        push!(preprocess, X, DataScaler) # ensure data is standardized before running PCA
        process!(X, preprocess.processors[end]) # process in place
        retval.X = copy(X)

        # Add a PCA scaler
        preprocess_pca = deepcopy(preprocess)
        n_PCA_features = min(nindicators, 20)
        push!(preprocess_pca, X, DataLinearTransform, n_PCA_features)

        Z = Array(Float64, n_PCA_features, nframes)
        retval.Z = process!(Z, X, preprocess_pca.processors[end])

        retval.Y = Y
        retval.extractor = extractor
        retval.preprocess = preprocess
        retval.preprocess_pca = preprocess_pca
        retval.action_clamper = DataClamper(
                Array(Float64, 2),
                vec(minimum(Y, 2)),
                vec(maximum(Y, 2))
            )

        retval
    end
end
type DynamicForestBehavior <: AbstractVehicleBehavior
    forest::Ensemble

    targets::ModelTargets
    indicators::Vector{AbstractVehicleBehavior}
    extractor::FeatureSubsetExtractor
    processor::DataPreprocessor
    action_clamper::DataClamper # contains action vector and clamping range

    # preallocated memory
    ϕ::Vector{Float64} # [n_predictors+1] (last one is always 1.0)

    function DynamicForestBehavior(
        forest::Ensemble,
        n_predictors::Int,
        params::DF_TrainParams,
        pre::DF_PreallocatedData,
        )

        extractor = pre.extractor
        action_clamper = pre.action_clamper
        processor = params.use_PCA ? pre.preprocess_pca : pre.preprocess

        retval = new()
        retval.forest = forest
        retval.targets = params.targets

        # copy the processor and ensure that the copied io arrays match up
        retval.extractor = FeatureSubsetExtractor(deepcopy(extractor.x), extractor.indicators)
        retval.processor = deepcopy(processor, retval.extractor)
        retval.action_clamper = DataClamper(Array(Float64, 2),
                                   deepcopy(action_clamper.f_lo),
                                   deepcopy(action_clamper.f_hi))

        retval.ϕ = Array(Float64, n_predictors+1)

        retval
    end
end
ProbDrive2016ITSC.get_targets(behavior::DynamicForestBehavior) = behavior.targets
ProbDrive2016ITSC.get_indicators(behavior::DynamicForestBehavior) = behavior.extractor.indicators

function preallocate_learning_data(
    dset::ModelTrainingData2,
    params::DF_TrainParams)

    DF_PreallocatedData(dset, params)
end

function _condition_predictor_mean!(
    normal::MvNormal,
    predictor_indeces::Vector{Int},
    A::Matrix{Float64},
    behavior::DynamicForestBehavior,
    )

    ϕ = behavior.ϕ

    for (i,p) in enumerate(predictor_indeces)
        ϕ[i] = behavior.processor.z[p]
    end
    ϕ[end] = 1.0

    normal.μ[:] = A * ϕ

    nothing
end
function select_action(
    behavior::DynamicForestBehavior,
    runlog::RunLog,
    sn::StreetNetwork,
    colset::UInt,
    frame::Int
    )

    observe!(behavior.extractor, runlog, sn, colset, frame)
    process!(behavior.processor)

    # pick a tree at random
    tree_index = rand(1:length(behavior.forest.trees))

    # sample from the MvNorm for said tree
    leaf = apply_tree(behavior.forest.trees[tree_index], behavior.processor.z)::AutoregressiveMvNormLeaf
    _condition_predictor_mean!(leaf.m, leaf.predictor_indeces, leaf.A, behavior)

    Distributions._rand!(leaf.m, behavior.action_clamper.x)
    process!(behavior.action_clamper)

    action_lat = behavior.action_clamper.x[1]
    action_lon = behavior.action_clamper.x[2]

    (action_lat, action_lon)
end

function _set_and_process_action!(behavior::DynamicForestBehavior, action_lat::Float64, action_lon::Float64)
    behavior.action_clamper.x[1] = action_lat
    behavior.action_clamper.x[2] = action_lon
    process!(behavior.action_clamper)
    behavior.action_clamper
end
function _calc_action_loglikelihood(behavior::DynamicForestBehavior)

    #=
    Action should be in behavior.action_clamper and should already have been processed
    Observations should be in behavior.processor.z and should already have been processed
    =#

    ntrees = length(behavior.forest.trees)
    log_ntrees = log(ntrees)

    total_probability_density = 0.0
    logl = 0.0
    for tree in behavior.forest.trees
        leaf = apply_tree(tree, behavior.processor.z)::AutoregressiveMvNormLeaf
        _condition_predictor_mean!(leaf.m, leaf.predictor_indeces, leaf.A, behavior)
        normal  = leaf.m

        p = pdf(normal, behavior.action_clamper.x)/ntrees

        if total_probability_density == 0.0
            logl = logpdf(normal, behavior.action_clamper.x) - log_ntrees
        elseif p < total_probability_density
            # NOTE: log(a + b) = log a + log(1 + b/a)
            logl = logl + log1p(p/total_probability_density)
        else
            # NOTE: log(a + b) = log b + log(1 + a/b)
            logl = (logpdf(normal, behavior.action_clamper.x) - log_ntrees) + log1p(total_probability_density/p)
        end
        total_probability_density += p
    end

    logl
end

function calc_action_loglikelihood(
    behavior::DynamicForestBehavior,
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
    behavior::DynamicForestBehavior,
    features::DataFrame,
    frameind::Integer,
    )

    observe!(behavior.extractor, features, frameind)
    process!(behavior.processor)

    action_lat = features[frameind, symbol(behavior.targets.lat)]::Float64
    action_lon = features[frameind, symbol(behavior.targets.lon)]::Float64
    _set_and_process_action!(behavior, action_lat, action_lon)

    _calc_action_loglikelihood(behavior)
end

function train(
    training_data::ModelTrainingData2,
    preallocated_data::DF_PreallocatedData,
    params::DF_TrainParams,
    foldset::FoldSet,
    )

    extractor = preallocated_data.extractor

    ntrees = params.ntrees
    max_tree_depth = params.max_tree_depth
    n_split_tries = params.n_split_tries
    min_samples_split = params.min_samples_split
    min_samples_leaves = params.min_samples_leaves
    n_random_predictor_samples = params.n_random_predictor_samples
    n_autoregression_predictors = params.n_autoregression_predictors
    min_split_improvement = params.min_split_improvement
    partial_sampling = params.partial_sampling
    autogression_coef = params.autogression_coef

    X_orig = params.use_PCA ? preallocated_data.Z : preallocated_data.X
    Y_orig = preallocated_data.Y

    # TODO(tim): do this without globals
    global DEFAULT_AUTOREGRESSION_CONSTANT = autogression_coef
    global DEFAULT_NUM_AUTOREGRESSION_PREDICTORS = n_autoregression_predictors
    global DEFAULT_NUM_PREDICTOR_SAMPLES = n_random_predictor_samples
    global GLOBAL_DIAG_SHRINKAGE_THRESOLD = params.diag_shrinkage_threshold

    build_tree_params = BuildTreeParameters(
        round(Int, sqrt(length(extractor.indicators))),
        max_tree_depth,
        min_samples_split,
        min_samples_leaves,
        min_split_improvement,
        n_split_tries,
        LossFunction_MSE,
        AutoregressiveMvNormLeaf
        )

    nframes = length(foldset)
    X = Array(Float64, nframes, size(X_orig, 1)) # NOTE: this is transpose of X_orig
    y = Array(Float64, 2, nframes)

    i = 0
    for frame in foldset
        i += 1

        for j in 1 : size(X_orig, 1)
            X[i,j] = X_orig[j,frame]
        end

        y[1,i] = Y_orig[1,frame]
        y[2,i] = Y_orig[2,frame]
    end
    @assert(i == nframes)

    ensemble = build_forest(y, X, ntrees, build_tree_params, partial_sampling)
    DynamicForestBehavior(ensemble, n_autoregression_predictors, params, preallocated_data)
end

end # end module
