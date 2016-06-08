module RandomForestBehaviors

using Reexport

include("MvDecisionTrees.jl")
using .MvDecisionTrees

include("GindeleRandomForestBehaviors.jl")
@reexport using .GindeleRandomForestBehaviors

include("DynamicForestBehaviors.jl")
@reexport using .DynamicForestBehaviors

end # module
