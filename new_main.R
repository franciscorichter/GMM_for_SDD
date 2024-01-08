#library(ggplot2)
#library(gridExtra)



N_max <- 100  # Number of species
true_alphas <- c(-1.5, -0.5, -0.01)  # Initial guess for alphas
true_betas <- c(0.05, 0.002, -0.03)  # Initial guess for betas
true_params <- c(true_alphas, true_betas)
max_t <- 100
max_attempts <- 50  # Maximum attempts for tree generation

result_obs <- create_tree_mat_fixed_N(N_max, true_alphas, true_betas, 0, max_attempts, max_t)

# Generate observed data for comparison
tree_mat_obs <- result_obs$tree_mat
tree_extant_obs <- L2phylo(unname(tree_mat_obs), dropextinct=TRUE)
stats_obs <- stats_PDM(tree_extant_obs)

plot(tree_extant_obs)

# Gradient Descent Settings
learn_rate <- 5e-1
patience <- 100
iters <- 1000
n_trees_D <- 10
n_trees_sgd <- 10


# Other necessary parameters
n_stats <- N_max * (N_max - 1) / 2  # Number of statistics
n_pars <- length(true_params)  # Number of parameters


# Function to generate randomized initial parameters for the search algorithm
generate_random_initial_params <- function(length, range) {
  # Generate random parameters within a given range
  return(runif(length, min = -range, max = range))
}

# Generate randomized initial parameters for the search algorithm
# The sign is not preserved and values can be both positive and negative
initial_alphas <- generate_random_initial_params(3, 1)  # Random alphas in range [-1, 1]
initial_betas <- generate_random_initial_params(3, 0.1) # Random betas in range [-0.1, 0.1]
initial_params <- c(initial_alphas, initial_betas)


# Run gradient descent
grad_descent_results <- gradDescentEvolutionModel(N_max, stats_obs, n_stats, n_pars, learn_rate, patience, 
                                                  iters, initial_params, n_trees_D, n_trees_sgd, max_t, max_attempts = max_attempts,print_inf=TRUE)


# Extract results
final_parameters <- grad_descent_results$pars[nrow(grad_descent_results$pars), ]
stats_diff_results <- grad_descent_results$stats_diff
times_taken <- grad_descent_results$times
# Iteration numbers
# Extracting parameter data from grad_descent_results
parameters <- grad_descent_results$pars

iterations <- 1:nrow(parameters)



# Plot time taken for each iteration
ggplot(data.frame(iteration=iterations[-c(1,2)], value=times_taken), aes(x=iteration, y=value)) +
  geom_line() + theme_bw() + 
  labs(title="Time Taken for Each Iteration", x="Iteration", y="Time (Seconds)")


# Number of parameters
n_params <- ncol(parameters)


# Create a list to store ggplot objects
plots <- list()

# Creating ggplot for each parameter
for (i in 1:n_params) {
  p <- ggplot(data.frame(iteration=iterations, value=parameters[, i]), aes(x=iteration, y=value)) +
    geom_line() + theme_bw() + 
    geom_hline(yintercept=true_params[i], color="red", linetype="dashed") +
    labs(title=paste("Parameter", i), x="Iteration", y="Value")
  plots[[i]] <- p
}

# Arrange the plots in a 2 by 3 grid
grid.arrange(grobs=plots, ncol=3, nrow=2)

ggplot(data.frame(iteration=iterations[-c(1,2)], value=grad_descent_results$stats_diff_abs), aes(x=iteration, y=value)) +
  geom_line() + theme_bw() + 
  labs(title="Mean Statsistics difference", x="Iteration", y="Value")


