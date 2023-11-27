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

# res contains the following:
  # pars: matrix with p columns corresponding to each parameter's value throughout the algorithm
  # stats_diff: matrix with columns being the individual statistics differences throghout the algorithm
  # stats_diff_abs: vector with the total absolute difference in statistics with observed data
  # times: the time for each iteration, which summed give the total runtime
  # D: a list of all the D matrices computed withing the algorithm, the size of which depends on how often we choose to calculate it


D_overlap = c()
N = c()
Ns = c(100000,200000)
for (N_trees in Ns){
  D1 = calc_D_stats(n_stats, n_pars, pars_i, N_trees, max_attempts, max_t, FALSE)
  D2 = calc_D_stats(n_stats, n_pars, pars_i, N_trees, max_attempts, max_t, FALSE)
  
  N = c(N, N_trees)
  ovr = sum((D1/D2)>0)/length(D1)
  D_overlap = c(D_overlap, ovr)
  print(N_trees)
}

plot(N, D_overlap, main="D overlap % vs number of trees", ylab="overlap")


stops = c(0,10,100,1000,5000,10000,20000,50000,100000)
#pars_i = true_pars
eps = pars_i/10

s = matrix(nrow=0,ncol=n_stats)
s1 = matrix(nrow=0, ncol=n_stats)
sp = matrix(nrow=0,ncol=n_stats)
s1p = matrix(nrow=0, ncol=n_stats)
D_test = list(list(),list(),list(),list(),list(),list(),list())


for (N in 1:length(stops)){
  
  # generate 1st original trees 
  for (i in (stops[N]+1):stops[N+1]){
    result = create_tree_mat_fixed_N(N_max, pars_i[1:3], pars_i[4:n_pars], 0, max_attempts, max_t)
    
    if (length(result)==0){
      return(c())
    }  
    
    if(result$attempt >= max_attempts){
      return(c())
    }
    
    tree_mat = result$tree_mat
    tree_extant = L2phylo(unname(tree_mat), dropextinct=TRUE)
    
    stats_obs_pdm = stats_PDM(tree_extant)
    
    s = rbind(s, stats_obs_pdm)
  }
  stats_orig = colMedians(s)
  
  # generate 2nd original trees 
  for (i in (stops[N]+1):stops[N+1]){
    result = create_tree_mat_fixed_N(N_max, pars_i[1:3], pars_i[4:n_pars], 0, max_attempts, max_t)
    
    if (length(result)==0){
      return(c())
    }  
    
    if(result$attempt >= max_attempts){
      return(c())
    }
    
    tree_mat = result$tree_mat
    tree_extant = L2phylo(unname(tree_mat), dropextinct=TRUE)
    
    stats_obs_pdm = stats_PDM(tree_extant)
    
    s1 = rbind(s1, stats_obs_pdm)
  }
  stats_orig1 = colMedians(s1)
  
  
  
  # generate 1st trees with tweaked pars for each par
  
  stats_grad = matrix(nrow=n_pars, ncol=n_stats)
  stats_grad1 = matrix(nrow=n_pars, ncol=n_stats)
  
  for (p in 1:n_pars){
    #print(p)
    pars_eps = pars_i
    pars_eps[p] = pars_eps[p] + eps[p]
    
    
    for (i in (stops[N]+1):stops[N+1]){
      result = create_tree_mat_fixed_N(N_max, pars_eps[1:3], pars_eps[4:n_pars], 0, max_attempts, max_t)
      
      if (length(result)==0){
        return(c())
      }  
      
      if(result$attempt >= max_attempts){
        return(c())
      }
      
      tree_mat = result$tree_mat
      tree_extant = L2phylo(unname(tree_mat), dropextinct=TRUE)
      
      stats_obs_pdm = stats_PDM(tree_extant)
      
      sp = rbind(sp, stats_obs_pdm)
    }
    # take average
    stats_avg = colMedians(sp)
    # populate row of matrix
    
    stats_grad[p,] = (stats_avg - stats_orig)/eps[p]
  }
  
  
  # generate 2nd trees with tweaked pars for each par
  for (p in 1:n_pars){
    #print(p)
    pars_eps = pars_i
    pars_eps[p] = pars_eps[p] + eps[p]
    
    
    for (i in (stops[N]+1):stops[N+1]){
      result = create_tree_mat_fixed_N(N_max, pars_eps[1:3], pars_eps[4:n_pars], 0, max_attempts, max_t)
      
      if (length(result)==0){
        return(c())
      }  
      
      if(result$attempt >= max_attempts){
        return(c())
      }
      
      tree_mat = result$tree_mat
      tree_extant = L2phylo(unname(tree_mat), dropextinct=TRUE)
      
      stats_obs_pdm = stats_PDM(tree_extant)
      
      s1p = rbind(s1p, stats_obs_pdm)
    }
    # take average
    stats_avg1 = colMedians(s1p)
    # populate row of matrix
    
    stats_grad1[p,] = (stats_avg1 - stats_orig1)/eps[p]
  }
  
  D_list = list(stats_grad, stats_grad1)
  D_test[[N]] = append(D_test[[N]], D_list)
  cat("Both Ds with ", N, " trees done\n", sep="")
}

# calculating overlap and plotting
D_overlap = c()
for (d in 1:6){
  ovr = sum((D_test[[d]][[1]]/D_test[[d]][[2]])>0)/36
  D_overlap = c(D_overlap, ovr)
}
plot(stops[2:7],D_overlap)


