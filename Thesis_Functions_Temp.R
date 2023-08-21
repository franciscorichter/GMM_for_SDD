library(DDD)
#library(treebalance)
#library(MASS)
library(ape)
#library(phytools)
#library(apTreeshape)
#library(castor)
#library(ggplot2)
#library(geiger)
library(numDeriv)
library(paleobuddy)

# calculate LTT statistic with number of timepoints = ltt_points
# tree is a phylo object
calc_ltt_stat = function(tree, ltt_points = 50){
  
  times = seq(round(tp-1),0,length.out=ltt_points)
  
  # calc ltt at all event times
  ltt_full = ape::ltt.plot.coords(tree)
  ltt_full[,"time"] = -1*ltt_full[,"time"]
  
  # calc ltt at arbitrary times
  ltt = c()
  for (t in times){# for arbitrary times
    c=1
    tt = ltt_full[c,"time"] # loop over all times
    while(t < tt){ # when t surpasses tt take the previous N value 
      c=c+1
      tt = ltt_full[c,"time"]
    }
    N = ltt_full[c-1,"N"]
    # bind to ltt  
    ltt = c(ltt,unname(N))
  }
  if(length(ltt)<length(times)){ltt=c(ltt,unname(ltt_full[c,"N"]))}
  return(ltt)
}


# initialization of the tree_mat
init_tree_mat = function(){
  #matrix with tree data
  tree_mat <- matrix(ncol = 4, nrow = 0)
  nomi <- c("birth.time", "parent", "node", "death.time")
  colnames(tree_mat) <- nomi
  # we start from 2 species 
  tree_mat = rbind(tree_mat, c(0, 0, 1, NA))
  tree_mat = rbind(tree_mat, c(0, 1, 2, NA))
  
  # list of lineages that are alive and not extinct (for event allocation)
  tips_at_t = c(1,2)
  # counter to add global lineage ids 
  counter = 2
  
  return(list(tree_mat = tree_mat, tips_at_t = tips_at_t, counter = 2))
}


########################### GRADIENT DESCENT ########################### 

# convert tree_mat to event matrix with (times,types)
tree_mat_to_events = function(tree_mat, tp){
  
  spec = cbind(unique(tree_mat[,1]), 1)
  ext = cbind(tree_mat[,4][tree_mat[,4]!=-1], 0)
  if(ncol(ext)<2){events = spec} else {events = rbind(spec,ext)}
  
  # MIGHT HAVE TO CONVERT
  events[,1] = tp - events[,1] 
  
  # put events in chrono order and loop through them
  events = events[order(events[,1], decreasing = FALSE),]
  #print("events created from tree_mat")
  return(events)
}


# Computes the W component (var-cov) of the D matrix for SGD
calc_W = function(n_stats, n_pars, stats, stats_avg, n_trees_D){
  
  res = matrix(data=0, nrow=n_stats, ncol=n_stats)
  z = matrix(data=0, nrow=n_stats, ncol=n_stats)
  
  # loop over simulated trees and find matrix for each, then average
  for (i in 1:n_trees_D){
    res = res + (stats[i,]-stats_avg)%*%t(stats[i,]-stats_avg)
  }
  
  res = (1/n_trees_D)*res
  
  # TAKE DIAGONAL
  d = diag(res)
  diag(z) = d
  
  # calc inverse
  W = solve(z)
  #W = ginv(res)
  #print("W computed")
  return(W)
}


# Computes the G component (likelihood gradients) of the D matrix for SGD
calc_G = function(n_stats, n_pars, stats, stats_obs, n_trees_D, likelihoods){
  # difference between simulated and observed stats - NB: TAKING DIFF BETWEEN MATRIX AND VECTOR, CHECK DIMS
  
  # TEST WITHOUT SUBTRACTING OBSERVED STATS FROM ALL ELEMENTS
  #stat_diff = sweep(stats,2,stats_obs) # dim = T x S
  stat_diff = stats
  
  # we need s x p matrix
  # transpose stat_diff to get S x T and multiply with T x P likelihoods
  # this should yield the S x P matrix required
  G = t(stat_diff)%*%likelihoods
  G = G/n_trees_D
  #print("G computed")
  return(G)
}


# Computes ideal matrix D
calc_ideal_D_new = function(tp, stats_obs, n_stats, n_pars, pars_i, n_trees_D, ltt_points, times, covariates_list, max_attempts){
  
  t_i = Sys.time()
  
  D = matrix(nrow=n_pars, ncol=n_stats)
  stats = matrix(nrow=0, ncol=n_stats)
  likelihoods = matrix(nrow=n_trees_D, ncol=n_pars)
  
  # simulate n_trees trees and compute their stats
  for (i in 1:n_trees_D){
    attempt = 0
    result = create_tree_mat_phy_COV(tp, pars_i[1], pars_i[2:length(pars_i)], covariates_list, attempt, max_attempts)
    tree_mat = result$tree_mat
    events = tree_mat_to_events(tree_mat, tp)

    # compute likelihoods
    lik <- calculate_gradients(pars_i, tree_mat, covariates_list, include_diversity)
    likelihoods[i,] = lik
    
    tree = L2phylo(unname(tree_mat), dropextinct=TRUE)
    # LTT STATS
    lt = calc_ltt_stat(tree, ltt_points)
    # ALL STATS
    s = c(lt)
    stats = rbind(stats, s)
    cat("D tree ", i, " created \n")
  }
  
  # round because LTT
  stats_avg = colMeans(stats)
  
  W = calc_W(n_stats, n_pars, stats, stats_avg, n_trees_D)
  G = calc_G(n_stats, n_pars, stats, stats_obs, n_trees_D, likelihoods)
  
  D = t(G)%*%W
  print("D computed")
  
  runtime_D = Sys.time() - t_i

  return(list(D=D, runtime_D=runtime_D))
}

# Calculates scaling matrix M
calc_M = function(D){
  DD = apply(D, 1, function(x) sd(x))
  M = diag(1/DD)
  return(M)
}

# GRADIENT DESCENT WITH CONSTANT D MATRIX
grad_descent_constD = function(tp, D, M, stats_obs, n_stats, n_pars, learn_rate, patience,
                               iters, pars_i, n_trees_D, n_trees_sgd, covariates_list){
  max_attempts=100
  # stats are stats of observed data
  c = 1 # iteration counter
  p = 0 # patience counter
  cond = Inf
  pars = matrix(nrow=0, ncol=length(pars_i))
  pars = rbind(pars, t(pars_i))
  pars_new = pars_i
  stats_diff = matrix(nrow=0,ncol=n_stats)
  stats_diff_abs = c()
  ll = learn_rate
  
  
  t_i = Sys.time()
  while (c < iters){ # CHANGE CONDITION

    stats_gen = matrix(nrow=0,ncol=n_stats)
    
    # n_trees_sgd: how many trees to generate at each step of sgd to calc stat diff
    for (i in 1:n_trees_sgd){
      attempt = 0
      result = create_tree_mat_phy_COV(tp, pars_new[1], pars_new[2:length(pars_new)], covariates_list,
                                       attempt, max_attempts)
      tree_mat = result$tree_mat
      tree = L2phylo(unname(tree_mat), dropextinct=TRUE)

      # LTT STATS
      lt = calc_ltt_stat(tree, ltt_points)
      # ALL STATS
      s = c(lt)
      stats_gen = rbind(stats_gen, s)
    }
    
    #print("end of sgd trees generation and stat calculation")
    stats_gen = colMeans(stats_gen)
    
    # compare stats
    diff = stats_gen - stats_obs 
    # update pars
    pars_new = pars_i - learn_rate*(M%*%D%*%(diff))
    
    # update condition (check diff in pars with prev iter and diff in stats on this iter)
    cond = pars_new - pars_i
    stats_diff = rbind(stats_diff, diff)
    stats_diff_abs = c(stats_diff_abs, sum(abs(diff)))
    
    # counter for LR reduction
    lr_count = 0
    # adaptive LR if B_N > 0 halve the LR and reiterate
    while (pars_new[length(pars_new)] > 0){
      #cat("LR: ", learn_rate, " -> ", learn_rate/2, "\n")
      learn_rate = learn_rate/2
      pars_new = pars_i - learn_rate*(M%*%D%*%(diff))
      lr_count = lr_count + 1
      
      # set limit on times we reduce lr_count
      if (lr_count > 100){
        cat("LR reduction limit reached: returning result \n")
        pars = rbind(pars, t(pars_new))
        runtime_SGD = Sys.time() - t_i
        learn_rate = ll
        return(list(pars=pars, stats_diff=stats_diff, stats_diff_abs=stats_diff_abs,
                    runtime_SGD=runtime_SGD, learn_rate=learn_rate))
      }
    }
    
    learn_rate = ll
    
    # if mu < 0 return/set to 0
    if (pars_new[1] < 0){
      pars_new[1] = 0
      # cat("mu < 0: returning result \n")
      # pars = rbind(pars, t(pars_new))
      # runtime_SGD = Sys.time() - t_i
      # return(list(pars=pars, stats_diff=stats_diff, stats_diff_abs=stats_diff_abs,
      #             runtime_SGD=runtime_SGD, learn_rate=learn_rate))
    }
    
    # initial lag to start checking for patience
    if (c>500){
      # check if stats_diff_abs[i] > stats_diff_abs[i-1] and if so increase patience
      target = mean(tail(stats_diff_abs, 20)[1:19]) # setting tail 20
      value = mean(tail(stats_diff_abs, 3)) # setting tail 3
      
      if (value > target){
        p = p + 1
        
        # check if p > max_p
        if (p > patience){
          cat("Reached patience limit of ", patience, ". Returning result \n")
          # remove last p values
          pars = head(pars, -p)
          stats_diff = head(stats_diff, -p)
          stats_diff_abs = head(stats_diff_abs, -p)
          # compute runtime
          runtime_SGD = Sys.time() - t_i
          # return result
          return(list(pars=pars, stats_diff=stats_diff, stats_diff_abs=stats_diff_abs,
                      runtime_SGD=runtime_SGD, learn_rate=learn_rate))
        }
      } else {p = 0}
    }
    
    
    # update pars_1 to pars_new
    pars_i = pars_new
    c = c + 1
    print(c)
    cat("pars_new = ", pars_i, "\n")
    cat("stats_diff = \n", round(diff,1), "\n")
    pars = rbind(pars, t(pars_i))
    
  }
  
  runtime_SGD = Sys.time() - t_i
  
  return(list(pars=pars, stats_diff=stats_diff, stats_diff_abs=stats_diff_abs,
              runtime_SGD=runtime_SGD, learn_rate=learn_rate))
}


# Runs SGD algorithm with arbitrary starting condition (default: mu=0.1, B0=0.75, BN=-0.002)
run_sgd_DD = function(tp, calc_D = TRUE, stats_obs, n_stats, n_pars, learn_rate, iters, patience,
                          pars_i, n_trees_D, n_trees_sgd, ltt_points, times, covariates_list){
  
  if (sum(calc_D) == TRUE){
    # calculating D matrix
    D_res = calc_ideal_D_new(tp, stats_obs, n_stats, n_pars, pars_i, n_trees_D, ltt_points, times, covariates_list, max_attempts)
    D = D_res$D
    runtime_D = round(D_res$runtime_D,2)
  } else {
    D = calc_D
    times = times
    runtime_D = 0
  }
  
  # calcuating M matrix
  M = calc_M(D)
  
  # performing sgd with IC = pars_i
  res = grad_descent_constD(tp, D, M, stats_obs, n_stats, n_pars, learn_rate, patience,
                            iters, pars_i, n_trees_D, n_trees_sgd, covariates_list)
  pars = res$pars
  stats_diff=res$stats_diff
  stats_diff_abs=res$stats_diff_abs
  runtime_SGD = round(res$runtime_SGD,2)
  runtime_tot = round(runtime_D+runtime_SGD,0)
  learn_rate = res$learn_rate
  
  return(list(pars=pars, stats_diff=stats_diff, stats_diff_abs=stats_diff_abs,
              D=D, M=M, runtime_tot=runtime_tot, learn_rate=learn_rate))
}



find_ic = function(tp, stats_obs, n_gen, n_treesXgen, n_stats, ltt_points, n_pars, pars_intervals, covariates_list,max_attempts){
  # matrices whcih will contain respective pars and mean(stats)
  stats = matrix(nrow=n_gen, ncol=n_stats)
  pars = matrix(nrow=n_gen, ncol=n_pars)
  
  for (i in 1:n_gen){
    
    # sample the pars for i-th generation form pars_intervals
    for (p in 1:n_pars){
      # pars_intervals should be (n_pars x 2) matrix with min and max
      pars[i,p] = runif(1, min=pars_intervals[p,1], max=pars_intervals[p,2])
    }
    
    
    # matrix that will contain stats from n_treesXgen trees in i-th generation
    stats_j = matrix(nrow=0, ncol=n_stats)
    # generate n_treesXgen trees, calc stats, take avg
    for(j in 1:n_treesXgen){
      # generate tree
      attempt=0
      res = create_tree_mat_phy_COV(tp, pars[i,1], pars[i,2:ncol(pars)], covariates_list, attempt, max_attempts)
      tree_mat = res$tree_mat
      # get extant tree
      tree_extant = DDD::L2phylo(unname(tree_mat), dropextinct=TRUE)
      # computing primary and LTT stats on observed tree
      stats_single = calc_ltt_stat(tree_extant, ltt_points)
      # append to stats
      stats_j = rbind(stats_j, stats_single)
    }
    
    # append mean to general stats
    stats[i,] = colMeans(stats_j)
    
    cat("find_ic(): gen ", i, "done\n", sep="")
  }
  
  # find abs() of difference in stats
  stat_diff = abs(sweep(stats, 2, stats_obs))
  # calc row sums
  sums = rowSums(stat_diff)
  # find minimum
  id_min = which(sums == min(sums))
  # set ic to pars[id_max,]
  ic = pars[id_min,]
    
  return(ic)
}



#################################################################################################


# function that computes the L matrix, which later will be converted to a phylo object
create_tree_mat_phy_COV = function(tp, mu_0, betas, covariates_list, attempt, max_attempts) {
  
  init = init_tree_mat()
  tree_mat = init$tree_mat
  tips_at_t = init$tips_at_t
  
  counter = init$counter
  # set t = 0, present time tp, total branch length
  t = 0
  left = c(1)
  right = c(2)
  
  # Gillespie algorithm to simulate birth-death process
  while (t < tp) {
    
    # allocate event to lineage by sampling from tips_at_t
    if (length(tips_at_t) != 0 & length(left) != 0 & length(right) != 0) {
      node = tips_at_t[sample.int(length(tips_at_t), 1)]
    } else { 
      cat("Generation interrupted as lineages have become extinct at time ", t, "\n")
      cat("Attempt n°", attempt, "\n")
      if (attempt >= max_attempts) {
        cat("Max number of attempts reached: ", max_attempts, "\n")
        stop()
      } else {
        
        sub_result <- create_tree_mat_phy_COV(tp, mu_0, betas, covariates_list, attempt+1, max_attempts)
        sub_tree_mat <- sub_result$tree_mat
        sub_attempt <- sub_result$attempt
        
        if (!is.null(sub_tree_mat)) {
          tree_mat <- sub_tree_mat
        }
        attempt <- sub_attempt
        
        return(list(tree_mat=tree_mat, attempt=attempt))
      }
    }
    
    # NEW VERSION
    #t_sample = get_nhpp_realization(lambda, tp, length(tips_at_t))
    
    # sample from NHPP
    t_sample = rexp.var(1, function(tm) length(tips_at_t)*(lambda(tm, covariates_list, length(tips_at_t), betas)+mu(mu_0)), 
             now=t, tMax=tp, shape=NULL, TS=0, fast=TRUE)
    
    # set new time
    t = t + t_sample
    #print(t)
    
    # break if over global timespan
    if (t > tp) {
      break
    }
    
    # Non-homogeneous exponential rate
    tot_rate = (mu(mu_0) + lambda(t, covariates_list, length(tips_at_t), betas))*length(tips_at_t)
    #cat("Total rate = ", tot_rate, "\n", sep="")
    
    
    if(tot_rate > 100000){
      cat("Tot_rate = " , tot_rate, "exceeded cap of 100000 \n")
      #cat("indiv spec rate = ", lambda_dd(t, tips_at_t, betas), "\n")
      #cat("indiv ext rate = ", mu_dd(mu), "\n")
      #cat("tot rate = ", tot_rate, "\n")
    }
    
    
    # sample event type
    sum_of_rates = mu(mu_0) + lambda(t, covariates_list, length(tips_at_t), betas)
    p = runif(1, 0, sum_of_rates)

    
    # if event speciation
    if (p > mu(mu_0)){
      # define childs and update counter
      child = counter+1
      counter = counter+1
      # insert node in corresponding lineage (left/right)
      if (node %in% left){left = c(left, child)}
      else if (node %in% right){right = c(right, child)}
      # insert edge in matrix and modify tips_at_t
      tree_mat = rbind(tree_mat, c(t, node, child, NA))
      tips_at_t = c(tips_at_t, child)
    } else {  # if extinction
      # remove node from corresponding lineage (left/right)
      if (node %in% left){left = left[left != node]}
      else if (node %in% right){right = right[right != node]}
      # interrupt the lineage by inserting death time and removing node from tips
      tree_mat[which(tree_mat[,3]==node), 4] = t
      tips_at_t = tips_at_t[tips_at_t != node]
    }
  }
  # put matrix in correct form by setting tp = 0 and ti = tp
  tree_mat[,"birth.time"] = tp - tree_mat[,"birth.time"]
  tree_mat[,"death.time"] = tp - tree_mat[,"death.time"]
  
  # set NAs = -1 for extant species
  tree_mat[is.na(tree_mat[,"death.time"]),"death.time"] = -1
  
  return(list(tree_mat=tree_mat, attempt=attempt))
}


# Define a function to evaluate covariate values
evaluate_covariates <- function(t, covariates) {
  closest_values <- sapply(covariates, function(covariate) {
    closest_time_index <- which.min(abs(covariate[,1] - t))
    return(covariate[closest_time_index, 2])
  })
  return(unname(closest_values))
}


# Define the lambda function using covariates
lambda <- function(t, covariates, N, betas) {
  if (length(covariates) == 0) {
    sp = max(0, betas[1] + betas[length(betas)] * N)
  } else {
    covariate_values <- evaluate_covariates(t, covariates)
    sp = max(0, betas[1] + sum(betas[2:(length(betas)-1)] * covariate_values) + betas[length(betas)] * N)
  }
  return(sp)
}

mu = function(mu_0){
  er = max(0, mu_0)
  return(er)
}

log_likelihood <- function(parameters, tree_data, covariates_list, include_diversity=TRUE) {
  # Extract parameters and betas from parameters
  mu_0 <- parameters[1]
  betas <- parameters[2:length(parameters)]
  #beta_diversity <- parameters[length(parameters)]
  
  log_lik <- 0
  events <- tree_mat_to_events(tree_data, tp)  # Convert tree data to events
  N <- 2  # Initialize diversity to 2
  
  for (e in 2:nrow(events)) {
    t = events[e, 1]
    t_prev = events[e - 1, 1]
    
    if (events[e, 2] == 1) {
      N <- N + 1  # Increase diversity for speciation
      #covariate_values <- evaluate_covariates(t, covariates_list)
      lambda_value <- lambda(t, covariates_list, N, betas)
      #print(lambda_value)  # Add this line
      log_lik = log_lik + log(lambda_value) - N * (lambda_value + mu(mu_0)) * (t - t_prev)
    } else {
      N <- N - 1  # Decrease diversity for extinction
      #covariate_values <- evaluate_covariates(t, covariates_list)
      lambda_value <- lambda(t, covariates_list, N, betas)
      #print(lambda_value)  # Add this line
      log_lik = log_lik + log(mu(mu_0)) - N * (lambda_value + mu(mu_0)) * (t - t_prev)
    }
    #cat("e = ", e, ", log-lik =", log_lik, "\n", sep="")
  }
  return(log_lik)
}


# Define function to calculate gradients
calculate_gradients <- function(parameters, tree_data, covariates_list, include_diversity) {
  gradient <- grad(func = log_likelihood, x = parameters, 
                   tree_data = tree_data, covariates_list = covariates_list, 
                   include_diversity = include_diversity)
  return(gradient)
}



# 
# tp=5
# mu_0=0.1
# betas=c(0.4, 0.05,-0.007)
# parameters = c(mu_0,betas)
# 
# n_stats=50
# n_pars=4
# learn_rate=1e-6
# iters = 10000
# pars_i = parameters
# n_trees_D = 100
# n_trees_sgd = 1
# ltt_points = 50
# times = seq(round(tp-1),0,length.out=ltt_points)
# max_attempts = 100
# patience=100
# 
# # Taking temperatures and simplifying
# library(paleobuddy)
# temp = paleobuddy::temp
# temp = temp[temp[,1]<=tp,]
# 
# temperatures = matrix(nrow=tp+1, ncol=2)
# temperatures[,1] = 0:tp
# for (r in 0:tp){
#   val = mean(temp[temp[,1] > (r-0.5) & temp[,1] < (r+0.5) ,2])
#   temperatures[r+1,2] = val
# }
# covs = list(temperatures = temperatures)
# 
# 
# # simulating extant tree
# result = create_tree_mat_phy_COV(tp, mu_0, betas, covs, 0, max_attempts)
# tree_mat = result$tree_mat
# tree_mat
# tree_extant = DDD::L2phylo(unname(tree_mat), dropextinct=TRUE)
# lt_obs = calc_ltt_stat(tree_extant, ltt_points)
# stats_obs = c(lt_obs)
# stats_obs
# 
# 
# D_res = calc_ideal_D_new(stats_obs,n_stats,n_pars,pars_i,n_trees_D,ltt_points,times,covs)
# D = D_res$D
# 
# ris = run_sgd_DD(calc_D = D, stats_obs, n_stats, n_pars, learn_rate, iters, patience,
#                       parameters, n_trees_D, n_trees_sgd, ltt_points, times, covs)
#   








