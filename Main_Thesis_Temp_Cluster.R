source("Thesis_Functions_Temp.R")


tp = 5 # total time (mln of years)
mu_0 = 0.1 # constant extinction rate
betas = c(0.4, 0.05, -0.007) # coefficients for speciation {B_0, B_1,...,B_T, B_N} (B_N should be last)
true_pars = c(mu_0, betas) 

n_stats = 50 # number of statistics (ltt points in this case)
n_pars = 4 # number of parameters (including mu)
learn_rate = 1e-5
iters = 3000 # max number of iterations for the code to stop
n_trees_D = 1000 # number of trees simulated to compute D
n_trees_sgd = 1 # number of trees simulated at each step of SGD
ltt_points = 50 # number of ltt points
times = seq(round(tp-1),0,length.out=ltt_points) # times
max_attempts = 100 # max attempts for full tree to be created (in case all species become extinct before the end, retry max_attempts times)
patience = 100 

# Intervals of pars to sample from
pars_intervals = matrix(data = c(c(0.01,0.3,0,-0.01), c(0.25,0.8,0.1,-0.003)), nrow=n_pars, ncol=2)
n_trees_ic = 500 # number of different ic combinations to try for find_ic()
n_trees_X_ic = 3 # number of trees to generate for each ic combination
  
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


# defining result list
ris_all = list()

  # populating tree_mat through Gillespie
  result = create_tree_mat_phy_COV(tp, mu_0, betas, covs, 0, max_attempts)
  tree_mat = result$tree_mat
  
  # converting to phylo and plotting
  #tree_full = DDD::L2phylo(unname(tree_mat), dropextinct=FALSE)
  tree_extant = DDD::L2phylo(unname(tree_mat), dropextinct=TRUE)
  
  # computing primary and LTT stats on observed tree
  stats_obs = calc_ltt_stat(tree_extant, ltt_points)
  # stats_obs = c(lt_obs)
  
  # find initial conditions
  pars_i = find_ic(tp, stats_obs, n_trees_ic, n_trees_X_ic, n_stats, ltt_points, n_pars, pars_intervals, covs, max_attempts)
  
  # true pars
  #pars_true = c(mu, betas)
  
  
  ris_all[[1]] = run_sgd_DD(tp, calc_D = TRUE, stats_obs, n_stats, n_pars, learn_rate, iters, patience,
                            pars_i, n_trees_D, n_trees_sgd, ltt_points, times, covs)
  

  
# # VIEW RESULT
# tail_size = 100
#   
# for (i in 1:n_pars){
#   val = round(mean(tail(ris_all[[1]]$pars[,i], tail_size)), 6)
#     
#   plot(1:nrow(ris_all[[1]]$pars), ris_all[[1]]$pars[,i], type="l",
#         main = paste("Par_", i, "=", val,".png"),
#         xlab = "Iters", ylab = paste("Par_", i))
#   }
  
  
# FOR SAVING RESULT: IF NOT ON CLUSTER MUST DEFINE J
cat("Tree n°", j, "simulated \n")
j <- as.integer(Sys.getenv("j"))
save.image(file = paste0("Sim_", j, ".RData"))




