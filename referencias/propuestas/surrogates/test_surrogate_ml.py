from surrogates_ml import *

# Set the ABM Evaluation Budget
budget = 500

# Set out-of-sample test and montecarlo sizes
test_size = 100
montecarlos = 100
calibration_threshold = 2.00

# Get an on out-of-sample test set that does not have combinations from the
# batch or iterative experiments
final_test_size = (test_size * montecarlos)

# Set the ABM parameters and support
islands_exploration_range = np.array([
        (0,10), # rho
        (0.8,2), # alpha
        (0.0,1.0), # phi
        (0.0,1.0), # pi
        (0.0,1.0), # eps
        (10,100), # N
        (0.0,1.0)]) # Lambda

n_dimensions = islands_exploration_range.shape[0]

# 1. Draw the Pool
# Set pool size
pool_size = 1000000

# Draw Pool
pool = get_sobol_samples(n_dimensions, pool_size, islands_exploration_range)

# 2. Randomly Draw and Evaluate Initialization Set
# Draw the initialization set from the pool as a permutation of the pool index.

# Set number of selections per round
samples_to_select = np.ceil(np.log(pool_size)).astype(int)

# Set initialization samples
initialization_samples = np.random.permutation(pool_size)[:samples_to_select]

# Evaluate Initialization Set
evaluated_set_X = pool[initialization_samples]
evaluated_set_y = evaluate_islands_on_set(evaluated_set_X)

# Update unevaluated_set_X
unevaluated_set_X = pool[list(set(range(pool_size)) - set(initialization_samples))]

surrogate_model, surrogate_parameter_space = set_surrogate_as_gbt()
print("Evaluated set size: ", evaluated_set_y.shape[0])

while evaluated_set_y.shape[0] < budget:
    print(evaluated_set_y.shape[0])
    # 3. Build Surrogate on evaluated samples
    surrogate_model_this_round = fit_surrogate_model( \
        evaluated_set_X,evaluated_set_y,
        surrogate_model=surrogate_model,
        surrogate_parameter_space=surrogate_parameter_space)

    # 4. Predict Response over Pool
    predict_response_pool = surrogate_model_this_round.predict(unevaluated_set_X)
    predicted_positives = calibration_condition(predict_response_pool,
                                                calibration_threshold)
    num_predicted_positives = predicted_positives.sum()

    # 5. Select small subset of Pool for Evaluation
    evaluated_set_X, evaluated_set_y, unevaluated_set_X = get_round_selections( \
            evaluated_set_X,evaluated_set_y,
            unevaluated_set_X,
            predicted_positives, num_predicted_positives,
            samples_to_select,
            budget)

# 6. Output Final Surrogate Model
surrogate_model = fit_surrogate_model(evaluated_set_X,evaluated_set_y,
                           surrogate_model=surrogate_model,
                           surrogate_parameter_space=surrogate_parameter_space)
