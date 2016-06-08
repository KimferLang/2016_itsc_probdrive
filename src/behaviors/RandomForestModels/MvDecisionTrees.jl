module MvDecisionTrees

export
    LossFunction,
    LossFunction_MSE,
    LossFunction_LOGL,
    LossFunction_LOGL_MEAN_SUBTRACTED,

    Leaf,
    MeanVecLeaf,
    CovLeaf,
    MvNormLeaf,

    InternalNode,
    Ensemble,
    BuildTreeParameters,

    depth,
    print_tree,
    apply_tree,
    apply_forest!,

    build_leaf,
    build_tree,
    build_forest,
    build_tree_parameters,

    get_split_index_set

import Distributions: MvNormal, fit_mle
import StatsBase: cov

if VERSION.minor ≥ 4
    typealias Range1{Int} Range{Int}
    _int(x) = round(Int, x, RoundNearestTiesUp)
else
    _int(x) = int(x)
end

##############################################################################
# TreeData

type TreeData{T<:AbstractFloat, U<:Real}
    # learning data
    Y::Matrix{T} # [n_observations × n_rows]
    X::Matrix{U} # [n_rows × n_predictors]
    assignment::Vector{Int} # [n_rows]

    # preallocated learning memory
    μ_l::Vector{T} # n
    μ_r::Vector{T} # n
    Σ_l::Matrix{T} # n×n
    Σ_r::Matrix{T} # n×n
    subfeature_indeces::Vector{Int}
end

Base.size(d::TreeData) = (size(d.Y,2), size(d.Y,1), size(d.X,2)) # (n, o, p)  n_rows, n_obs, n_predictors
get_num_samples(d::TreeData) = size(d.Y,2)
get_num_observations(d::TreeData) = size(d.Y,1)
function get_id_count(data::TreeData, id::Int)
    id_count = 0
    for i = 1 : size(data.X, 1)
        if data.assignment[i] == id
            id_count += 1
        end
    end
    id_count
end
function get_nth_valid_id(data::TreeData, id::Int, n::Int)
    id_count = 0
    for i = 1 : size(data.X, 1)
        if data.assignment[i] == id
            id_count += 1
            if id_count == n
                return i
            end
        end
    end
    return 0
end
function get_left_right_counts{T,U}(data::TreeData{T,U}, predictor_index::Int, id::Int, threshold::T)

    num_left = 0
    num_right = 0

    for i = 1 : get_num_samples(data)
        if data.assignment[i] == id
            if data.X[i,predictor_index] < threshold
                num_left += 1
            else
                num_right += 1
            end
        end
    end

    (num_left, num_right)
end

##############################################################################
# LossFunction

abstract LossFunction
type LossFunction_MSE <: LossFunction end # mean squared error loss function
type LossFunction_LOGL <: LossFunction end # log likelihood loss function
type LossFunction_LOGL_MEAN_SUBTRACTED <: LossFunction end # log likelihood loss function where labels already has predicted mean subtracted

function loss{T<:AbstractFloat, U<:Real}(::Type{LossFunction_MSE},
    data::TreeData{T,U},
    predictor_index::Int,
    assignment_id::Int,
    threshold::T, # split point
    )

    # mse loss for a vector is tr{E{(y-p)(y-p)ᵀ}}
    # where y is the label
    # and p is the prediction - the average of labels on that side of the threshold

    Σ_l = data.Σ_l
    Σ_r = data.Σ_r

    n, o, p = size(data) # n_rows, n_obs, n_predictors
    mse = 0.0
    nl, nr = get_left_right_counts(data, predictor_index, assignment_id, threshold)

    for i = 1 : n
        if data.assignment[i] == assignment_id
            s_l, s_r = 0.0, 0.0
            for j = 1 : o
                l = data.Y[j,i]
                mse += l*l
                if data.X[j, predictor_index] < threshold
                    s_l += l
                else
                    s_r += l
                end
            end
            mse -= s_l*s_l / nl
            mse -= s_r*s_r / nr
        end
    end

    return (-mse, nl, nr)
end
function loss{T<:AbstractFloat, U<:Real}(::Type{LossFunction_MSE},
    data::TreeData{T,U},
    assignment_id::Int,
    )

    # mse loss for a vector is tr{E{(y-p)(y-p)ᵀ}}
    # where y is the label
    # and p is the prediction - the average of labels on that side of the threshold

    Σ = data.Σ_l
    n, o, p = size(data) # n_samples, n_outcomes, n_predictors
    mse = 0.0

    id_count = get_id_count(data, assignment_id)

    for i = 1 : n
        if data.assignment[i] == assignment_id
            s = 0.0
            for j = 1 : o
                l = data.Y[j,i]
                mse += l*l
                s += l
            end
            mse -= s*s / id_count
        end
    end

    return -mse
end

function loss{T<:AbstractFloat, S<:Real}(::Type{LossFunction_LOGL},
    data::TreeData{T,S},
    predictor_index::Int,
    assignment_id::Int,
    threshold::T, # split point
    )

    μ_l = data.μ_l
    μ_r = data.μ_r
    Σ_l = data.Σ_l
    Σ_r = data.Σ_r

    n, o, p = size(data)
    nl, nr = 0, 0

    # compute the mean
    for i = 1 : o
        μ_l[i] = 0.0
        μ_r[i] = 0.0
    end
    for i = 1 : n
        if data.assignment[i] == assignment_id
            if data.X[i,predictor_index] < threshold
                nl += 1
                for a = 1 : o
                    μ_l[a] += data.Y[a,i]
                end
            else
                nr += 1
                for a = 1 : o
                    μ_r[a] += data.Y[a,i]
                end
            end
        end
    end
    if nl < 2 || nr < 2
        # NOTE(tim): force at least two samples on each side to ensure proper covariance matrix
        return (-Inf, nl, nr)
    end
    for i = 1 : o
        μ_l[i] /= nl
        μ_r[i] /= nr
    end

    # solve for upper triangle of symmetric covariance matrix
    for i = 1 : o*o
        Σ_l[i] = 0.0
        Σ_r[i] = 0.0
    end
    for i = 1 : n
        if data.assignment[i] == assignment_id
            if data.X[i,predictor_index] < threshold
                for a = 1 : o
                    l = data.Y[a,i] - μ_l[a]
                    for b = a : o
                        Σ_l[a,b] += l*(data.Y[b,i] - μ_l[b])
                    end
                end
            else
                for a = 1 : o
                    l = data.Y[a,i] - μ_r[a]
                    for b = a : o
                        Σ_r[a,b] += l*(data.Y[b,i] - μ_r[b])
                    end
                end
            end
        end
    end
    for i = 1 : o*o
        Σ_l[i] /= (nl-1)
        Σ_r[i] /= (nr-1)
    end
    for a = 2:o
        # copy over symmetric component
        for b = 1:a-1
            Σ_l[a,b] = Σ_l[b,a]
            Σ_r[a,b] = Σ_r[b,a]
        end
    end

    logdetΣ_l = log(max(_fast_det(Σ_l), 1e-20))
    logdetΣ_r = log(max(_fast_det(Σ_r), 1e-20))

    # compute log likelihood score under gaussian distribution
    # NOTE: does not compute the true logl likelihood
    #       computation ignores constant components
    logl = 0.0
    for i = 1 : n
        if data.assignment[i] == assignment_id
            if data.X[i] < threshold
                logl -= logdetΣ_l
                for a = 1 : o
                    subval = 0.0
                    for b = 1 : o
                        subval += Σ_l[a,b]*(data.Y[b,i] - μ_l[b])
                    end
                    logl -= subval*(data.Y[a,i] - μ_l[a])
                end
            else
                logl -= logdetΣ_r
                for a = 1 : o
                    subval = 0.0
                    for b = 1 : o
                        subval += Σ_r[a,b]*(data.Y[b,i] - μ_r[b])
                    end
                    logl -= subval*(data.Y[a,i] - μ_r[a])
                end
            end
        end
    end

    return (logl, nl, nr)
end
function loss{T<:AbstractFloat, S<:Real}(::Type{LossFunction_LOGL},
    data::TreeData{T,S},
    assignment_id::Int,
    )

    μ = data.μ_l
    Σ = data.Σ_l

    n, o, p = size(data)
    # n, m = size(data.Y)
    m_assigned = 0

    # compute the sample mean
    for i = 1 : o
        μ[i] = 0.0
    end
    for i = 1 : n
        if data.assignment[i] == assignment_id
            m_assigned += 1
            for a = 1 : o
                μ[a] += data.Y[a,i]
            end
        end
    end
    if m_assigned ≤ 2
        # insufficient entries to split
        return -Inf
    end
    for i = 1 : o
        μ[i] /= m_assigned
    end

    # compute covariance matrix

    for i = 1 : o*o
        Σ[i] = 0.0
    end
    for i = 1 : n
        if data.assignment[i] == assignment_id
            for a = 1 : o
                l = data.Y[a,i] - μ[a]
                for b = a : o
                    Σ[a,b] += l*(data.Y[b,i] - μ[b])
                end
            end
        end
    end
    for i = 1 : o*o
        Σ[i] /= (m_assigned-1)
    end
    for a = 2:o
        # copy over symmetric component
        for b = 1:a-1
            Σ[a,b] = Σ[b,a]
        end
    end

    logdetΣ = log(max(_fast_det(Σ), 1e-20))

    # compute log likelihood score under gaussian distribution
    # NOTE: does not compute the true logl likelihood
    #       computation ignores constant components
    logl = 0.0
    for i = 1 : n
        if data.assignment[i] == assignment_id
            logl -= logdetΣ
            for a = 1 : o
                subval = 0.0
                for b = 1 : o
                    subval += Σ[a,b]*(data.Y[b,i] - μ[b])
                end
                logl -= subval*(data.Y[a,i] - μ[a])
            end
        end
    end
    logl
end

function loss{T<:AbstractFloat, S<:Real}(::Type{LossFunction_LOGL_MEAN_SUBTRACTED},
    data::TreeData{T,S},
    predictor_index::Int,
    assignment_id::Int,
    threshold::T, # split point
    )

    # data.Y is n_samples × n_outcomes and contains the distance from the mean, x-μ
    # we want to predict the covariance matrix Σ, which is determined by the covariance
    # of the samples in a leaf

    Σ_l = data.Σ_l
    Σ_r = data.Σ_r

    n, o, p = size(data)

    nl, nr = 0, 0

    # solve for upper triangle of symmetric matrix
    for i = 1 : o*o
        Σ_l[i] = 0.0
        Σ_r[i] = 0.0
    end

    for i = 1 : n
        if data.assignment[i] == assignment_id
            if data.X[i,predictor_index] < threshold
                for a = 1 : o
                    l = data.Y[a,i]
                    for b = a : o
                        Σ_l[a,b] += l*data.Y[b,i]
                    end
                end
                nl += 1
            else
                for a = 1 : o
                    l = data.Y[a,i]
                    for b = a : o
                        Σ_r[a,b] += l*data.Y[b,i]
                    end
                end
                nr += 1
            end
        end
    end
    for i = 1 : o*o
        Σ_l[i] /= (nl-1)
        Σ_r[i] /= (nr-1)
    end

    # NOTE(tim): force at least two samples on each side to ensure proper covariance matrix
    if nl < 2 || nr < 2
        return (-Inf, nl, nr)
    end

    # copy over symmetric component
    for a = 2:o
        for b = 1:a-1
            Σ_l[a,b] = Σ_l[b,a]
            Σ_r[a,b] = Σ_r[b,a]
        end
    end

    logdetΣ_l = log(max(_fast_det(Σ_l), 1e-20))
    logdetΣ_r = log(max(_fast_det(Σ_r), 1e-20))

    # compute log likelihood (we seek to maximize it)
    # NOTE: ignores constant components
    logl = 0.0
    for i = 1 : n
        if data.assignment[i] == assignment_id
            if data.X[i] < threshold
                logl -= logdetΣ_l
                for a = 1 : o
                    subval = 0.0
                    for b = 1 : o
                        subval += Σ_l[a,b]*data.Y[b,i]
                    end
                    logl -= subval*data.Y[a,i]
                end
            else
                logl -= logdetΣ_r
                for a = 1 : o
                    subval = 0.0
                    for b = 1 : o
                        subval += Σ_r[a,b]*data.Y[b,i]
                    end
                    logl -= subval*data.Y[a,i]
                end
            end
        end
    end

    return (logl, nl, nr)
end
function loss{T<:AbstractFloat, S<:Real}(::Type{LossFunction_LOGL_MEAN_SUBTRACTED},
    data::TreeData{T,S},
    assignment_id::Int,
    )

    # data.Y is n_samples × n_outcomes and contains the distance from the mean, x-μ
    # we want to predict the covariance matrix Σ, which is determined by the covariance
    # of the samples in a leaf

    Σ = data.Σ_l
    n, o, p = size(data)

    for i = 1 : o*o
        Σ[i] = 0
    end

    # solve for upper triangle of symmetric matrix

    m_assigned = 0
    for i = 1 : n
        if data.assignment[i] == assignment_id
            m_assigned += 1
            for a = 1 : o
                l = data.Y[a,i]
                for b = a : o
                    Σ[a,b] += l*data.Y[b,i]
                end
            end
        end
    end
    if m_assigned ≤ 2
        return -Inf
    end

    for i = 1 : o*o
        Σ[i] /= (m_assigned-1)
    end

    # copy over symmetric component
    for a = 2:o
        for b = 1:a-1
            Σ[a,b] = Σ[b,a]
        end
    end

    logdetΣ = log(max(_fast_det(Σ), 1e-20))

    # compute log likelihood score under gaussian distribution
    # NOTE: does not compute the true logl likelihood,
    #       computation ignores constant components
    logl = 0.0
    for i = 1 : n
        if data.assignment[i] == assignment_id
            logl -= logdetΣ
            for a = 1 : o
                subval = 0.0
                for b = 1 : o
                    subval += Σ[a,b]*data.Y[b,i]
                end
                logl -= subval*data.Y[a,i]
            end
        end
    end
    logl
end

##############################################################################
# Leaf

abstract Leaf
immutable MeanVecLeaf{T<:AbstractFloat} <: Leaf
    μ::Vector{T}
end
immutable CovLeaf{T<:AbstractFloat} <: Leaf
    Σ::Matrix{T}
end
immutable MvNormLeaf <: Leaf
    m::MvNormal
end

function build_leaf{T<:AbstractFloat, S<:Real}(::Type{MeanVecLeaf},
    data::TreeData{T,S},
    assignment_id::Int,
    )

    #=
    Construct a leaf whose payload is the mean vector of samples within it
    =#

    n, o, p = size(data) # n_rows, n_obs, n_predictors
    leaf_mean = zeros(T, o)

    n_for_id = 0
    for i = 1 : n
        if data.assignment[i] == assignment_id
            n_for_id += 1
            for j = 1 : o
                leaf_mean[j] += data.Y[j,i]
            end
        end
    end
    for j = 1 : o
        leaf_mean[j] /= n_for_id
    end

    MeanVecLeaf(leaf_mean)
end
function build_leaf{T<:AbstractFloat, S<:Real}(::Type{CovLeaf},
    data::TreeData{T,S},
    assignment_id::Int,
    )

    #=
    Construct a leaf whose payload is the covariance matrix for all samples within it
    =#
    o = get_num_observations(data)
    CovLeaf(_calc_covariance!(zeros(T, o, o), data, assignment_id))
end
function build_leaf{T<:AbstractFloat, S<:Real}(::Type{MvNormLeaf},
    data::TreeData{T,S},
    assignment_id::Int,
    )

    #=
    Construct a leaf whose payload is the mvnormal for all samples within it
    =#

    n, o, p = size(data)
    leaf_mean = zeros(T, o)

    n_for_id = 0
    for i = 1 : n
        if data.assignment[i] == assignment_id
            n_for_id += 1
            for j = 1 : o
                leaf_mean[j] += data.Y[j,i]
            end
        end
    end
    for j = 1 : o
        leaf_mean[j] /= n_for_id
    end

    # calc covariance

    Σ = zeros(T, o, o)

    for i = 1 : n
        if data.assignment[i] == assignment_id
            for a = 1 : o
                l = data.Y[a,i] - leaf_mean[a]
                for b = a : o
                    Σ[a,b] += l*(data.Y[b,i]-leaf_mean[b])
                end
            end
        end
    end


    if n_for_id < 2
        warn("MvDecisionTrees.build_leaf{MvNormLeaf}: Not enough samples ($n_for_id) to build a covariance matrix")
        for i = 1 : o
            Σ[i,i] = 1.0
        end
    else
        for i = 1 : o*o
            Σ[i] /= (n_for_id-1)
        end
        _diagonal_shrinkage!(Σ)
        for a = 2:o
            for b = 1: a-1
                Σ[a,b] = Σ[b,a]
            end
        end
    end

    for a = 2:o
        for b = 1:a-1
            Σ[a,b] = Σ[b,a]
        end
    end

    MvNormLeaf(MvNormal(leaf_mean, Σ))
end

##############################################################################

immutable InternalNode{T<:Any, LeafType<:Leaf}
    featid::Int # id of the feature we use to split
    featval::T # value over which we split
    left::Union{LeafType, InternalNode}
    right::Union{LeafType, InternalNode}
end

immutable Ensemble{LeafType<:Leaf}
    trees::Vector{Union{InternalNode, LeafType}}
end

type BuildTreeParameters{F<:LossFunction, LeafType<:Any}
    nsubfeatures::Int # number of candidate features to select for split; 0 will use all
    max_depth::Int # will only split if above this depth (depth of just a leaf is 0)
    min_samples_split::Int # minimum number of samples requried in a node to split
    min_samples_leaves::Int # minimum number of samples in newly created leaves
    min_split_improvement::Float64 # will only split if loss(split(leaf)) > loss(leaf) + min_split_improvement
    n_split_tries::Int # number of random splits to try (sampled uniformally from available predictors and data points)
    loss_function::Type{F}
    leaf_type::Type{LeafType} # the type used for the leaf
end
function build_tree_parameters{F<:LossFunction, LeafType<:Leaf}(;
    nsubfeatures::Int=0,
    max_depth::Int=typemax(Int),
    min_samples_split::Int=5,
    min_samples_leaves::Int=0,
    min_split_improvement::Float64=-Inf,
    loss_function::Type{F}=LossFunction_MSE,
    leaf_type::Type{LeafType}=MeanVecLeaf,
    )

    @assert(nsubfeatures ≥ 0)
    @assert(max_depth ≥ 0)
    @assert(min_samples_split ≥ 0)
    @assert(min_samples_leaves ≤ min_samples_split/2)

    BuildTreeParameters{F, LeafType}(nsubfeatures, max_depth, min_samples_split,
               min_samples_leaves, min_split_improvement, loss_function, leaf_type)
end

function Base.print(io::IO, params::BuildTreeParameters)
    println(io, "BuildTreeParameters:")
    @printf(io, "\tnsubfeatures:          %8d\n", params.nsubfeatures)
    @printf(io, "\tmax_depth:             %8d\n", params.max_depth)
    @printf(io, "\tmin_samples_split:     %8d\n", params.min_samples_split)
    @printf(io, "\tmin_samples_leaves:    %8d\n", params.min_samples_leaves)
    @printf(io, "\tmin_split_improvement: %8.3f\n", params.min_split_improvement)
    @printf(io, "\tloss_function:         %8s\n", string(params.loss_function))
end

Base.length(::Leaf) = 1
Base.length(tree::InternalNode) = length(tree.left) + length(tree.right)
Base.length(ensemble::Ensemble) = length(ensemble.trees)

depth(::Leaf) = 0
depth(tree::InternalNode) = 1 + max(depth(tree.left), depth(tree.right))

function print_tree(leaf::Leaf, depth::Integer=-1, indent::Integer=0)
    # NOTE(tim): I think this only works for classification
    print("    " ^ indent * "R-> ")
    println(leaf)
end
function print_tree(tree::InternalNode, depth::Integer=-1, indent::Integer=0)
    if depth == indent
        println()
        return
    end
    println("Feature $(tree.featid), Threshold $(tree.featval)")
    print("    " ^ indent * "L-> ")
    print_tree(tree.left, depth, indent + 1)
    print("    " ^ indent * "R-> ")
    print_tree(tree.right, depth, indent + 1)
end

apply_tree(leaf::Leaf, feature::Vector) = leaf
function apply_tree(tree::InternalNode, features::Vector)
    if features[tree.featid] < tree.featval
        apply_tree(tree.left, features)
    else
        apply_tree(tree.right, features)
    end
end

function apply_forest!{T<:Any}(μ::Vector{T}, forest::Ensemble{MeanVecLeaf}, features::Vector)

    # return the mean among tree's observations

    ntrees = length(forest)
    o = length(μ)

    for i in 1 : o
        μ[i] = 0.0
    end
    for i in 1 : ntrees
        leaf = apply_tree(forest.trees[i], features)::MeanVecLeaf
        for j in 1 : o
            μ[j] += leaf.μ[j]
        end
    end
    for i in 1 : o
        μ[i] /= ntrees
    end
    μ
end
function apply_forest!{T<:Any}(Σ::Matrix{T}, forest::Ensemble{CovLeaf}, features::Vector)

    # return the mean covariance matrix among the tree's covariances

    ntrees = length(forest)
    len = length(Σ)

    for i in 1 : len
        Σ[i] = 0.0
    end
    for i in 1 : ntrees
        leaf = apply_tree(forest.trees[i], features)::CovLeaf
        for j in 1 : len
            Σ[j] += leaf.Σ[j]
        end
    end
    for i in 1 : len
        Σ[i] /= ntrees
    end

    Σ
end

function _split{F<:LossFunction, T<:AbstractFloat, U<:Real}(
    data::TreeData{T,U},
    assignment_id::Int,

    min_samples_leaves::Int,
    min_split_improvement::Float64,
    loss_function::Type{F},
    n_thresholds::Int,
    )

    # returns a tuple: (index_of_feature_we_split_over, threshold)

    nrow, nfeatures = size(data.X)
    best = (0,0)
    best_loss = loss(loss_function, data, assignment_id) + min_split_improvement
    n_thresholds = min(nrow-1, n_thresholds)

    _reservoir_sample!(data.subfeature_indeces, nfeatures)

    for i in data.subfeature_indeces

        f_lo = Inf
        f_hi = -Inf
        for j in 1 : nrow
            if data.assignment[j] == assignment_id
                v = data.X[j,i]
                if v < f_lo
                    f_lo = v
                elseif v > f_hi
                    f_hi = v
                end
            end
        end

        for j in 1 : n_thresholds

            threshold = f_lo + (f_hi - f_lo)*(j/(n_thresholds+1))

            value, num_left, num_right = loss(loss_function, data, i, assignment_id, threshold)

            if value > best_loss && num_left ≥ min_samples_leaves && num_right ≥ min_samples_leaves
                best_loss = value
                best = (i, threshold)
            end
        end
    end

    best
end
function _split_unif_rand{F<:LossFunction, T<:AbstractFloat, U<:Real}(
    data::TreeData{T,U},
    assignment_id::Int,

    min_samples_leaves::Int,
    min_split_improvement::Float64,
    loss_function::Type{F},
    n_split_tries::Int,
    )

    # returns a tuple: (index_of_feature_we_split_over, threshold)

    id_count = get_id_count(data, assignment_id)
    nrows, nfeatures = size(data.X)
    id_range = 1:id_count
    predictor_range = 1:nfeatures

    @assert(id_count > 0)
    @assert(nfeatures > 0)

    best_id = 0
    best_threshold = 0.0
    best_loss = loss(loss_function, data, assignment_id) + min_split_improvement

    for _ in 1 : n_split_tries
        predictor_id = rand(predictor_range)
        nth_id = rand(id_range)
        row_index = get_nth_valid_id(data, assignment_id, nth_id)
        threshold = data.X[row_index, predictor_id]
        value, num_left, num_right = loss(loss_function, data, predictor_id, assignment_id, threshold)

        if value > best_loss && num_left ≥ min_samples_leaves && num_right ≥ min_samples_leaves
            best_loss = value
            best_id = predictor_id
            best_threshold = threshold
        end
    end

    (best_id, best_threshold)
end

function _build_tree{T<:AbstractFloat, U<:Real}(
    data::TreeData{T,U},
    assignment_id::Int,
    params::BuildTreeParameters,
    _depth::Int=0,
    )

    n = get_num_samples(data)
    n_for_id = 0
    for i = 1 : n
        n_for_id += convert(Int, data.assignment[i] == assignment_id)
    end

    # println("size(assignment): ", size(assignment))
    # println(assignment[1:20])
    # println("assignment_id:    ", assignment_id)
    # println("n_for_id: ", n_for_id)
    # println("_depth ≥ params.max_depth   ", _depth ≥ params.max_depth)
    # println("n_for_id < params.min_samples_split:   ", n_for_id < params.min_samples_split)

    if _depth ≥ params.max_depth || n_for_id < params.min_samples_split
        return build_leaf(params.leaf_type, data, assignment_id)
    end

    id, threshold = _split(data, assignment_id,
                        params.min_samples_leaves,
                        params.min_split_improvement,
                        params.loss_function,
                        params.n_split_tries,
                        )

    # println("split id:      ", id)
    # println("split thresh:  ", threshold)

    if id == 0
        return build_leaf(params.leaf_type, data, assignment_id)
    end


    next_id = data.assignment[end]+1
    data.assignment[end] = next_id

    for i in 1 : n
        if data.assignment[i] == assignment_id && data.X[i,id] < threshold
            data.assignment[i] = next_id
        end
    end

    return InternalNode{T, params.leaf_type}(id, threshold,
                _build_tree(data, next_id,       params, _depth+1),
                _build_tree(data, assignment_id, params, _depth+1))
end
function build_tree{T<:AbstractFloat, U<:Real}(labels::Matrix{T}, features::Matrix{U}, params::BuildTreeParameters)

    o, n = size(labels)
    μ_l = Array(T, o)
    μ_r = Array(T, o)
    Σ_l = Array(T, o,o)
    Σ_r = Array(T, o,o)

    nfeatures = size(features,2)
    nsubfeatures = params.nsubfeatures > 0 ? min(params.nsubfeatures, nfeatures) : nfeatures
    subfeature_indeces = Array(Int, nsubfeatures)
    assignment = ones(Int, n+1)

    data = TreeData{T,U}(labels, features, assignment, μ_l, μ_r, Σ_l, Σ_r, subfeature_indeces)

    _build_tree(data, 1, params)
end

function _build_forest_parallel{T<:AbstractFloat, U<:Real}(labels::Matrix{T}, features::Matrix{U}, ntrees::Integer, params::BuildTreeParameters, partialsampling::Float64)
    nlabels = _nsamples_labels(labels)
    nsamples = _int(partialsampling * nlabels)

    forest = @parallel (vcat) for i in 1:ntrees
        inds = rand(1:nlabels, nsamples)
        build_tree(labels[:,inds], features[inds,:], params)
    end
    Ensemble{params.leaf_type}(convert(Vector{Union{InternalNode, params.leaf_type}}, forest))
end
function _build_forest_serial{T<:AbstractFloat, U<:Real}(labels::Matrix{T}, features::Matrix{U}, ntrees::Integer, params::BuildTreeParameters, partialsampling::Float64)

    o, nlabels = size(labels)
    μ_l = Array(T, o)
    μ_r = Array(T, o)
    Σ_l = Array(T, o,o)
    Σ_r = Array(T, o,o)

    nfeatures = size(features,2)
    nsubfeatures = params.nsubfeatures > 0 ? min(params.nsubfeatures, nfeatures) : nfeatures
    subfeature_indeces = Array(Int, nsubfeatures)
    assignment = ones(Int, nlabels+1)
    nsamples = _int(partialsampling * nlabels)

    chosen_samples = Array(Int, nsamples)

    data = TreeData{T,U}(labels, features, assignment, μ_l, μ_r, Σ_l, Σ_r, subfeature_indeces)

    forest = Array(Union{InternalNode, params.leaf_type}, ntrees)
    for i = 1 : ntrees

        # select nsamples
        _reservoir_sample!(chosen_samples, nlabels)
        for j = 1:nlabels
            assignment[j] = 0 # remove labels that we are not considering
        end
        for j in chosen_samples
            assignment[j] = 1
        end
        assignment[end] = 1


        # build the tree
        forest[i] = _build_tree(data, 1, params)
    end

    Ensemble{params.leaf_type}(forest)
end
function build_forest{T<:AbstractFloat, U<:Real}(labels::Matrix{T}, features::Matrix{U}, ntrees::Integer, params::BuildTreeParameters, partialsampling::Float64=0.7)

    @assert(0.0 < partialsampling ≤ 1.0)
    if nprocs() == 1
        _build_forest_serial(labels, features, ntrees, params, partialsampling)
    else
        _build_forest_parallel(labels, features, ntrees, params, partialsampling)
    end
end

#####################################
# analysis functions

get_split_index_set(leaf::Leaf) = Set{Int}()
function get_split_index_set(tree::InternalNode)
    # returns a Set{Int} of indeces used as splits in the tree
    retval = Set{Int}()
    push!(retval, tree.featid)
    for featid in get_split_index_set(tree.left)
        push!(retval, featid)
    end
    for featid in get_split_index_set(tree.right)
        push!(retval, featid)
    end
    retval
end
function get_split_index_set{L<:Leaf}(forest::Ensemble{L})
    # returns a Set{Int} of indeces used as splits in the forest
    retval = Set{Int}()
    for tree in forest.trees
        for featid in get_split_index_set(tree)
            push!(retval, featid)
        end
    end
    retval
end

#####################################
# helper functions

function _calc_covariance!{T, U}(
    Σ::Matrix{T},
    data::TreeData{T,U},
    assignment_id::Int,
    )

    # println(fit_mle(MvNormal, data.Y))

    n, o, p = size(data)
    leaf_mean = zeros(T, o)

    n_for_id = 0
    for i = 1 : n
        if data.assignment[i] == assignment_id
            n_for_id += 1
            for j = 1 : o
                leaf_mean[j] += data.Y[j,i]
            end
        end
    end
    for j = 1 : o
        leaf_mean[j] /= n_for_id
    end

    # calc covariance

    Σ = zeros(T, o, o)

    for i = 1 : n
        if data.assignment[i] == assignment_id
            for a = 1 : o
                l = data.Y[a,i] - leaf_mean[a]
                for b = a : o
                    Σ[a,b] += l*(data.Y[b,i]-leaf_mean[b])
                end
            end
        end
    end

    for i = 1 : o*o
        Σ[i] /= (n_for_id-1) # we assume that we estimated the mean, so we use (n-1)
    end

    _diagonal_shrinkage!(Σ)

    # copy over symmetric part of matrix
    for a = 2:o
        for b = 1:a-1
            Σ[a,b] = Σ[b,a]
        end
    end

    Σ
end
function _diagonal_shrinkage!{T<:Real}(
    Σ::Matrix{T}, # 2×2 matrix with upper diagonal filled in
                  # this method ONLY affects the upper diagonal
    ε::T = convert(T, 1e-6),
    )

    # now we want to use Shrinkage Estimation - https://en.wikipedia.org/wiki/Estimation_of_covariance_matrices#Shrinkage_estimation
    # Since with a small number of samples we might get a Σ that isn't positive semidefinite
    # A matrix is pos. semidefinite it all of its eigenvalues are > 0
    # We choose γ such that A = γΣ + (1-γI) with all of eigenvalues are > ε
    # This method will work in general, but what is coded below assumes Σ is 2×2
    #      Σ = [a b; b d]
    @assert(size(Σ,1)==2 && size(Σ,2)==2)

    a, b, d = Σ[1,1], Σ[1,2], Σ[2,2]

    γ = 2.0*(1.0-ε) / (2.0 - a - d + sqrt(a*a -2*a*d + 4*b*b + d*d))

    if 0.0 < γ < 1.0
        # @assert(γ ≥ 0.0)

        Σ[1,1] = γ*a + (1.0-γ)
        Σ[1,2] = γ*b
        Σ[2,2] = γ*d + (1.0-γ)
    end

    Σ
end
function _fast_det{T<:AbstractFloat}(M::AbstractMatrix{T})
    n = size(M,1)
    if n == 1
        M[1]
    elseif n == 2
        M[1]*M[4] - M[2]*M[3]
    elseif n == 3
        M[1]*M[5]*M[9] + M[4]*M[8]*M[3] + M[7]*M[2]*M[6] -
            M[7]*M[5]*M[3] - M[4]*M[2]*M[9] - M[1]*M[8]*M[6]
    else
        det(M)
    end
end
function _reservoir_sample!(indeces::Vector{Int}, n::Int, k::Int=length(indeces))

    # performs Reservoir Sampling to obtain k elements from [1:n]
    #  indeces = input array
    #  n       = length of set of indeces to choose from (ie, we select from [1:n])
    #  k       = number of samples to choose and place in first k elements of indeces
    # runtime: O(n)
    # https://en.wikipedia.org/wiki/Reservoir_sampling

    if n == 0
        n = k
    end

    for i = 1 : k
        indeces[i] = i
    end

    for i = k+1 : n
        j = rand(1:i)
        if j ≤ k
            indeces[j] = i
        end
    end

    indeces
end

end # end module
