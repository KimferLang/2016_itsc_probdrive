module GindeleRandomForestBehaviors

using ProbDrive2016ITSC

using RandomForestBehaviors.MvDecisionTrees

import ProbDrive2016ITSC:
    AbstractVehicleBehavior,
    AbstractVehicleBehaviorPreallocatedData,
    AbstractVehicleBehaviorTrainParams,

    preallocate_learning_data,
    select_action,
    calc_action_loglikelihood,
    train

export
    GindeleRandomForestBehavior,
    GRF_TrainParams,

    select_action,
    calc_action_loglikelihood,
    train

type GRF_TrainParams <: AbstractVehicleBehaviorTrainParams

    targets::ModelTargets
    indicators::Vector{AbstractFeature}

    ntrees::Int # number of trees in the forest
    max_tree_depth::Int # max depth for trees in the forest
    min_samples_split::Int # min num samples in a node for split to be considered
    min_samples_leaves::Int # min num samples in resulting leaves for split to be done
    n_split_tries::Int # number of random split attempts
    use_PCA::Bool

    min_split_improvement::Float64 # minimum improvement for split to be accepted
    partial_sampling::Float64 # percentage of dataset to sample

    function GRF_TrainParams(;
        targets::ModelTargets = ModelTargets(Features.FUTUREDESIREDANGLE, Features.FUTUREACCELERATION),
        indicators::Vector{AbstractFeature} = AbstractFeature[],

        ntrees::Integer=10,
        max_tree_depth::Integer=5,
        min_samples_split::Integer=100,
        min_samples_leaves::Integer=20,
        n_split_tries::Integer=10,
        use_PCA::Bool = false,
        min_split_improvement::Float64=0.0,
        partial_sampling::Float64=0.7,
        )

        retval = new()
        retval.targets = targets
        retval.indicators = indicators
        retval.ntrees = ntrees
        retval.max_tree_depth = max_tree_depth
        retval.min_samples_split = min_samples_split
        retval.min_samples_leaves = min_samples_leaves
        retval.n_split_tries = n_split_tries
        retval.use_PCA = use_PCA
        retval.min_split_improvement = min_split_improvement
        retval.partial_sampling = partial_sampling
        retval
    end
end
type GRF_PreallocatedData <: AbstractVehicleBehaviorPreallocatedData

    # the full dataset matrices pre-extracted
    X::Matrix{Float64} # [nfeatures × nframes]
    Y::Matrix{Float64} # [2 × nframes]
    Z::Matrix{Float64} # [n_pca_features × nframes]

    extractor::FeatureSubsetExtractor
    preprocess::DataPreprocessor
    preprocess_pca::DataPreprocessor
    action_clamper::DataClamper

    function GRF_PreallocatedData(dset::ModelTrainingData2, params::GRF_TrainParams)
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

        extractor = FeatureSubsetExtractor(Array(Float64, nindicators), indicators)
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
type GindeleRandomForestBehavior <: AbstractVehicleBehavior
    model_μ::Ensemble
    model_Σ::Ensemble

    targets::ModelTargets
    indicators::Vector{AbstractFeature}
    extractor::FeatureSubsetExtractor
    processor::DataPreprocessor
    action_clamper::DataClamper

    # preallocated memory
    μ::Vector{Float64} # [2]
    Σ::Matrix{Float64} # [2×2]

    function GindeleRandomForestBehavior(
        model_μ::Ensemble,
        model_Σ::Ensemble,
        params::GRF_TrainParams,
        pre::GRF_PreallocatedData,
        )

        extractor = pre.extractor
        action_clamper = pre.action_clamper
        processor = params.use_PCA > 0 ? pre.preprocess_pca : pre.preprocess

        retval = new()
        retval.model_μ = model_μ
        retval.model_Σ = model_Σ

        retval.targets = params.targets
        # copy the processor and ensure that the copied io arrays match up
        retval.extractor = FeatureSubsetExtractor(deepcopy(extractor.x), extractor.indicators)
        retval.processor = deepcopy(processor, retval.extractor)
        retval.action_clamper = DataClamper(Array(Float64, 2),
                                   deepcopy(action_clamper.f_lo),
                                   deepcopy(action_clamper.f_hi))

        retval.μ = Array(Float64, 2)
        retval.Σ = Array(Float64, 2, 2)

        retval
    end
end
ProbDrive2016ITSC.get_targets(behavior::GindeleRandomForestBehavior) = behavior.targets
function ProbDrive2016ITSC.get_indicators(behavior::GindeleRandomForestBehavior)

    # run through and pull the indicators that are actually used in the tree

    indicators = Set{AbstractFeature}()

    for index in get_split_index_set(behavior.model_μ)
        push!(indicators, behavior.extractor.indicators[index])
    end
    for index in get_split_index_set(behavior.model_Σ)
        push!(indicators, behavior.extractor.indicators[index])
    end

    return collect(indicators)
end

function preallocate_learning_data(
    dset::ModelTrainingData2,
    params::GRF_TrainParams)

    GRF_PreallocatedData(dset, params)
end

function _calc_mvnormal(behavior::GindeleRandomForestBehavior)

    apply_forest!(behavior.μ, behavior.model_μ, behavior.processor.z)
    apply_forest!(behavior.Σ, behavior.model_Σ, behavior.processor.z)

    MvNormal(behavior.μ, behavior.Σ)
end

function select_action(
    behavior::GindeleRandomForestBehavior,
    runlog::RunLog,
    sn::StreetNetwork,
    colset::UInt,
    frame::Int
    )

    observe!(behavior.extractor, runlog, sn, colset, frame)
    process!(behavior.processor)

    normal = _calc_mvnormal(behavior)

    Distributions._rand!(normal, behavior.action_clamper.x)
    process!(behavior.action_clamper)

    action_lat = behavior.action_clamper.x[1]
    action_lon = behavior.action_clamper.x[2]

    (action_lat, action_lon)
end

function calc_action_loglikelihood(
    behavior::GindeleRandomForestBehavior,
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

    behavior.action_clamper.x[1] = action_lat
    behavior.action_clamper.x[2] = action_lon

    observe!(behavior.extractor, runlog, sn, colset, frame)
    process!(behavior.processor)

    normal = _calc_mvnormal(behavior)
    logpdf(normal, behavior.action_clamper.x)
end
function calc_action_loglikelihood(
    behavior::GindeleRandomForestBehavior,
    features::DataFrame,
    frameind::Integer,
    )

    behavior.action_clamper.x[1] = features[frameind, symbol(behavior.targets.lat)]::Float64
    behavior.action_clamper.x[2] = features[frameind, symbol(behavior.targets.lon)]::Float64

    observe!(behavior.extractor, features, frameind)
    process!(behavior.processor)

    apply_forest!(behavior.μ, behavior.model_μ, behavior.processor.z)
    apply_forest!(behavior.Σ, behavior.model_Σ, behavior.processor.z)

    mvnorm_model = MvNormal(behavior.μ, behavior.Σ)
    logpdf(mvnorm_model, behavior.action_clamper.x)
end

function train(
    training_data::ModelTrainingData2,
    preallocated_data::GRF_PreallocatedData,
    params::GRF_TrainParams,
    foldset::FoldSet,
    )

    ntrees = params.ntrees
    max_tree_depth = params.max_tree_depth
    min_samples_split = params.min_samples_split
    min_samples_leaves = params.min_samples_leaves
    n_split_tries = params.n_split_tries
    min_split_improvement = params.min_split_improvement
    partial_sampling = params.partial_sampling

    #######

    build_tree_params = BuildTreeParameters(
        ntrees, max_tree_depth, min_samples_split, min_samples_leaves,
        min_split_improvement, n_split_tries, LossFunction_MSE, MeanVecLeaf)

    X_orig = params.use_PCA ? preallocated_data.Z : preallocated_data.X

    Y_orig = preallocated_data.Y

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

    μ = build_forest(y, X, ntrees, build_tree_params, partial_sampling)

    mean_vector = Array(Float64, 2)
    for i = 1 : size(y, 1)
        y[:,i] -= apply_forest!(mean_vector, μ, vec(X[i,:])) # subtract the predicted mean
    end

    build_tree_params = BuildTreeParameters(
        ntrees, max_tree_depth, min_samples_split, min_samples_leaves,
        min_split_improvement, n_split_tries, LossFunction_MSE, CovLeaf)

    Σ = build_forest(y, X, ntrees, build_tree_params, partial_sampling)

    GindeleRandomForestBehavior(μ, Σ, params, preallocated_data)
end

end # end module
