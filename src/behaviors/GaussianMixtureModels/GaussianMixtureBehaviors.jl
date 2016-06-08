module GaussianMixtureBehaviors

using Reexport

include("GaussianMixtureRegressors.jl")
@reexport using .GaussianMixtureRegressors

end # module
