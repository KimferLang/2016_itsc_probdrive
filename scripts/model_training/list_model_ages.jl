using ProbDrive2016ITSC
using RandomForestBehaviors
using DynamicBayesianNetworkBehaviors

##############################
# PARAMETERS
##############################

include(Pkg.dir("ProbDrive2016ITSC", "scripts", "extract_params.jl"))

##############################
# PARAMETERS
##############################

context_classes = ["freeflow", "following", "lanechange"]
ncontext_classes = length(context_classes)

model_names = AbstractString[]
push!(model_names, "Static Gaussian")
push!(model_names, "Linear Gaussian")
push!(model_names, "Random Forest")
push!(model_names, "Dynamic Forest")
push!(model_names, "Mixture Regression")
push!(model_names, "Bayesian Network")
push!(model_names, "Linear Bayesian")
push!(model_names, "Static Gaussian Clean")
push!(model_names, "Linear Gaussian Clean")
push!(model_names, "Random Forest Clean")
push!(model_names, "Dynamic Forest Clean")
push!(model_names, "Mixture Regression Clean")
push!(model_names, "Bayesian Network Clean")
push!(model_names, "Linear Bayesian Clean")

current_time = now()

@printf("%-30s %-15s %4s %5s\n", "model name", "context class", "fold", "days")

for model_name in model_names

    println(uppercase(model_name))
    model_output_name = replace(lowercase(model_name), " ", "_")    # ex: bayesian_network
    model_short_name = convert_model_name_to_short_name(model_name) # ex: BN

    for context_class in context_classes

        dset_filepath_modifier = "_" * context_class

        for file in readdir(EVALUATION_DIR)
            if startswith(file, model_short_name * dset_filepath_modifier * "_fold")

                fold = parse(Int, match(r"(?<=_fold)\d+", file).match)

                model_for_fold_path_jld = joinpath(EVALUATION_DIR, file)
                model_save_time = JLD.load(model_for_fold_path_jld, "time")
                elapsed_days = round(Int, (current_time - model_save_time).value / 86400000)

                @printf("%-30s %-15s %2d     %3d\n", model_name, context_class, fold, elapsed_days)
            end
        end
    end

    println("")
end

println("DONE")