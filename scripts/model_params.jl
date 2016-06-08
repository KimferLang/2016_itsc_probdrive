behaviorset = Dict{AbstractString, BehaviorTrainDefinition}()

###################
# Partial Parameter Set

behaviorset["Static Gaussian"] = BehaviorTrainDefinition(SG_TrainParams())
behaviorset["Linear Gaussian"] = BehaviorTrainDefinition(
                                            LG_TrainParams(indicators=INDICATOR_SET),
                                            [
                                                BehaviorParameter(:ridge_regression_constant, [0.0,1e-5,1e-4,1e-3], 2)
                                            ])
behaviorset["Random Forest"] = BehaviorTrainDefinition(
                                            GRF_TrainParams(indicators=INDICATOR_SET),
                                            [
                                                BehaviorParameter(:ntrees, [3,10], 2),
                                                BehaviorParameter(:max_tree_depth, [2,3,4], 2),
                                                BehaviorParameter(:use_PCA, [false, true], 1),
                                                BehaviorParameter(:n_split_tries, [2,10,20,30], 2),
                                                BehaviorParameter(:min_split_improvement, [0.0,1.0], 1),
                                                BehaviorParameter(:partial_sampling, [0.2,0.5,0.6,0.8], 2),
                                            ])
behaviorset["Dynamic Forest"] = BehaviorTrainDefinition(
                                            DF_TrainParams(indicators=INDICATOR_SET),
                                            [
                                                BehaviorParameter(:ntrees, [3,10,20,30], 2),
                                                BehaviorParameter(:max_tree_depth, [2,3,4], 2),
                                                BehaviorParameter(:use_PCA, [false, true], 1),
                                                BehaviorParameter(:n_split_tries, [10,20,30], 2),
                                                BehaviorParameter(:min_split_improvement, [0.0,1.0], 2),
                                                BehaviorParameter(:partial_sampling, [0.2,0.5,0.6,0.8], 2),
                                            ])
behaviorset["Bayesian Network"] = BehaviorTrainDefinition(
                                            BN_TrainParams(
                                                indicators=INDICATOR_SET,
                                                preoptimize_target_bins=true,
                                                preoptimize_indicator_bins=true,
                                                optimize_structure=true,
                                                optimize_target_bins=false,
                                                optimize_parent_bins=false,
                                                ncandidate_bins=20,
                                                max_parents=5,
                                                nbins_lat=20,
                                                nbins_lon=10,
                                                dirichlet_prior=BDeuPrior(0.5),
                                            ),
                                            [
                                                BehaviorParameter(:max_parents, [3], 1),
                                                BehaviorParameter(:nbins_lat, 6:8, 2),
                                                BehaviorParameter(:nbins_lon, 6:8, 2),
                                            ]
                                            )
behaviorset["Linear Bayesian"] = BehaviorTrainDefinition(
                                            LB_TrainParams(indicators=INDICATOR_SET),
                                            [
                                                BehaviorParameter(:ridge_regression_constant, [0.0, 0.01, 0.1, 0.25], 2),
                                                BehaviorParameter(:min_σ_lat, [1e-7, 1e-6, 1e-5], 2),
                                                BehaviorParameter(:min_σ_lon, [1e-7, 5e-6, 1e-6], 2),
                                                BehaviorParameter(:max_parents, [2,3,4,5], 1),
                                            ])
behaviorset["Mixture Regression"] = BehaviorTrainDefinition(
                                            GMR_TrainParams(indicators=INDICATOR_SET, n_gmm_iter=20, n_init=1, min_covar=1e-7, n_components=5, max_n_indicators=3, unlearned_component_weight=0.05),
                                            BehaviorParameter[
                                                BehaviorParameter(:n_components, [3,4,5], 2),
                                                BehaviorParameter(:tol, [0.1, 0.15, 0.2], 2),
                                                BehaviorParameter(:min_covar, [1e-7,1e-6,1e-5,1e-4], 2),
                                                BehaviorParameter(:max_n_indicators, [3,4,5], 2),
                                                BehaviorParameter(:unlearned_component_weight, [0.01,0.05,0.1], 2)
                                                BehaviorParameter(:use_PCA, [false, true], 1),
                                            ])