source("Functions_N_MAX_lam&mu(N,P).R")

N_max = 50 # number of species
max_t = 100 # maximum time
alphas = c(-2, -1, -0.03) # coefficients for extinction {A_0, A_P, A_N}
betas = c(0.1, 0.001, -0.05) # coefficients for speciation {B_0, B_P, B_N}
true_pars = c(alphas, betas) 



n_stats = N_max*(N_max-1)/2 # number of statistics
n_pars = length(true_pars) # number of parameters
learn_rate = 1e-1 # learning rate
iters = 3000 # max number of iterations for SGD to stop
n_trees_D = 500 # number of trees simulated to compute D
steps_D = 500 # D is calculated every steps_D iterations
n_trees_sgd = 5 # number of trees simulated at each step of SGD
max_attempts = 50 # max attempts for full tree to be created (in case all species become extinct before the end, retry max_attempts times)
patience = 100 # stopping criterion (NOT IMPLEMENTED)


# Simulating the data
result = create_tree_mat_fixed_N(N_max, alphas, betas, 0, 100, max_t)
tree_mat = result$tree_mat

# converting to phylo and plotting
tree_full = DDD::L2phylo(unname(tree_mat), dropextinct=FALSE)
plot(tree_full)
tree_extant = DDD::L2phylo(unname(tree_mat), dropextinct=TRUE)
plot(tree_extant)

# Compute the statistics (phylogenetic distance matrix)
stats_obs = stats_PDM(tree_extant)

# Setting initial conditions
pars_i = true_pars

# Running SGD with MoM
res = grad_descent_constD(N_max, stats_obs, n_stats, n_pars, learn_rate, patience,
                          iters, pars_i, n_trees_D, n_trees_sgd, max_t, steps_D, print_info=TRUE)




