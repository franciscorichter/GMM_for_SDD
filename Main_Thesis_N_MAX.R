source("Thesis_Functions_N_MAX.R")

N_max = 100
tp = 5 # total time (mln of years)
mu_0 = log(0.1) # constant extinction rate
betas = c(0.5, 0.1, -0.01) # coefficients for speciation {B_0, B_1,...,B_T, B_N} (B_N should be last)
true_pars = c(mu_0, betas) 

n_stats = N_max*(N_max-1)/2 # number of statistics (ltt points in this case)
n_pars = 4 # number of parameters (including mu)
learn_rate = 1e-6
iters = 5000 # max number of iterations for the code to stop
n_trees_D = 500 # number of trees simulated to compute D
n_trees_sgd = 1 # number of trees simulated at each step of SGD
ltt_points = 100 # number of ltt points
#times = seq(round(tp-0.1),0,length.out=ltt_points) # times
max_attempts = 100 # max attempts for full tree to be created (in case all species become extinct before the end, retry max_attempts times)
patience = 100 

# Intervals of pars to sample from
pars_intervals = matrix(data = c(c(0.01,0.3,0,-0.01), c(0.25,0.8,0.1,-0.003)), nrow=n_pars, ncol=2)
n_trees_ic = 500 # number of different ic combinations to try for find_ic()
n_trees_X_ic = 2 # number of trees to generate for each ic combination

# Taking temperatures and simplifying
temp = paleobuddy::temp
temp = temp[temp[,1]<=tp,]

temperatures = matrix(nrow=tp+1, ncol=2)
temperatures[,1] = 0:tp
for (r in 0:tp){
  val = mean(temp[temp[,1] > (r-0.5) & temp[,1] < (r+0.5) ,2])
  temperatures[r+1,2] = val
}

# defining covariates
covs = list(temperatures = temperatures)


# populating tree_mat through Gillespie
result = create_tree_mat_fixed_N(N_max, mu_0, betas, covs, 0, 100)
tree_mat = result$tree_mat
tree_mat
# extant tree mat
tree_mat_extant = tree_mat[which(tree_mat[,"death.time"]==-1),]
tree_mat_extant
# total time
t_tot = result$t_tot
t_tot

# converting to phylo and plotting
tree_full = DDD::L2phylo(unname(tree_mat), dropextinct=FALSE)
plot(tree_full)
tree_extant = DDD::L2phylo(unname(tree_mat), dropextinct=TRUE)
plot(tree_extant)


# Compute the phylogenetic distance matrix
stats_obs = stats_PDM(tree_extant)

# computing primary and LTT stats on observed tree
# stats_obs = calc_ltt_stat(tree_extant, ltt_points)
# stats_obs
# stats_obs = c(lt_obs)




# Create an example tree
tree_newick <- "(((A:1,B:1):1,C:2):1,(D:1.5,(E:0.5,F:0.5):1):1.5);"
tree_extant <- read.tree(text = tree_newick)
plot(tree_extant)

PDM=calc_PDM(tree_extant)
PDM_r=reorder_PDM(PDM)
vectorize_PDM(PDM_r)


