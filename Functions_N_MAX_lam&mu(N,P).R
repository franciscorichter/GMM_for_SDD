library(DDD)
library(ape)
library(numDeriv)
library(paleobuddy)
library(tidyverse)
library(matrixStats)


# calculate LTT statistic with number of timepoints = ltt_points
# tree is a phylo object
calc_ltt_stat = function(tree, ltt_points = 50){
  
  
  # calc ltt at all event times
  ltt_full = ape::ltt.plot.coords(tree)
  ltt_full[,"time"] = -1*ltt_full[,"time"]
  
  tp = ltt_full[1,1]
  times = seq(tp-(tp/20),0,length.out=ltt_points)
  
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


# Calculates scaling matrix M
calc_M = function(D){
  DD = apply(D, 1, function(x) sd(x))
  M = diag(1/DD)
  return(M)
}

# GRADIENT DESCENT WITH CONSTANT D MATRIX
grad_descent_constD = function(N_max, stats_obs, n_stats, n_pars, learn_rate, patience,
                               iters, pars_i, n_trees_D, n_trees_sgd, max_t, steps_D, print_info=TRUE){
  max_attempts=50
  # stats are stats of observed data
  c = 0 # iteration counter
  p = 0 # patience counter
  cond = Inf
  pars = matrix(nrow=0, ncol=length(pars_i))
  pars = rbind(pars, t(pars_i))
  pars_new = pars_i
  stats_diff = matrix(nrow=0,ncol=n_stats)
  stats_diff_abs = c()
  times = c()
  #MD = D
  #MD = sweep(MD, 1, abs(exp(pars_i)), '*')
  
  # Set scaler (parameters at this point)
  scaling = abs(pars_new)
  
  D_info = list()
  
  
  while (c < iters){ # CONDITION
    t_i = Sys.time()  
    
    if (c %% steps_D == 0){
      D = calc_D_stats(n_stats,n_pars,pars_i,n_trees_D,max_attempts,max_t)
      
      if(length(D)==0){D = D_info[[length(D_info)]]}
      
      #M = calc_M(D)
      #MD = M%*%D
      MD = D/rowSums(abs(D))
      #MD = calc_M(D)%*%D
      #MD = MD/rowSums(abs(MD))
      
      # Store D matrix info
      D_info = rbind(D_info, list(D))
      
    }
    
    stats_gen = matrix(nrow=0,ncol=n_stats)
    
    # n_trees_sgd: how many trees to generate at each step of sgd to calc stat diff
    for (i in 1:n_trees_sgd){
      attempt = 0
      result = create_tree_mat_fixed_N(N_max, pars_new[1:3], pars_new[4:length(pars_new)],
                                          attempt, max_attempts, max_t)
      if (length(result)==0){
        return(list(pars=pars, stats_diff=stats_diff, stats_diff_abs=stats_diff_abs, D=D, times=times))
      }  
      
      att = result$attempt
      
      # if max_attempts is reached fix increment and change
      mu_incr = 0.1*abs(pars_new[1:3])
      betas_incr = 0.1*abs(pars_new[4:length(pars_new)])
      
      while (att >= max_attempts){
        pars_new[1:3] = pars_new[1:3] - mu_incr
        pars_new[4:length(pars_new)] = pars_new[4:length(pars_new)] + betas_incr
        result = create_tree_mat_fixed_N(N_max, pars_new[1:3], pars_new[4:length(pars_new)],
                                         attempt, max_attempts, max_t)
        if (length(result)==0){
          return(list(pars=pars, stats_diff=stats_diff, stats_diff_abs=stats_diff_abs, D=D, times=times))
        }
        att = result$attempt
        pars_i = pars_new
      }
      
      tree_mat = result$tree_mat
      tree = L2phylo(unname(tree_mat), dropextinct=TRUE)
      
      # PDM & LTT STATS
      s_pdm = stats_PDM(tree)
      stats_gen = rbind(stats_gen, s_pdm)
    }
    
    #print("end of sgd trees generation and stat calculation")
    stats_gen = colMedians(stats_gen)
    
    # compare stats
    diff = stats_gen - stats_obs  
    
    # WEIGHTS: USED TO RENORMALIZE DIFF (NB ALSO * N_STATS AFTER)
    #weights = abs(diff)/sum(abs(diff))
    
    # update pars
    pars_new = pars_i - learn_rate*(MD%*%(as.matrix(diff)))*scaling
    
    # update condition (check diff in pars with prev iter and diff in stats on this iter)
    cond = pars_new - pars_i
    stats_diff = rbind(stats_diff, diff)
    stats_diff_abs = c(stats_diff_abs, sum(abs(diff)))
    
    
    # initial lag to start checking for patience
    # if (c>500){
    #   # check if stats_diff_abs[i] > stats_diff_abs[i-1] and if so increase patience
    #   target = mean(tail(stats_diff_abs, 20)[1:19]) # setting tail 20
    #   value = mean(tail(stats_diff_abs, 3)) # setting tail 3
    #   
    #   if (value > target){
    #     p = p + 1
    #     
    #     # check if p > max_p
    #     if (p > patience){
    #       cat("Reached patience limit of ", patience, ". Returning result \n")
    #       # remove last p values
    #       pars = head(pars, -p)
    #       stats_diff = head(stats_diff, -p)
    #       stats_diff_abs = head(stats_diff_abs, -p)
    #       # compute runtime
    #       runtime_SGD = Sys.time() - t_i
    #       # return result
    #       return(list(pars=pars, stats_diff=stats_diff, stats_diff_abs=stats_diff_abs,
    #                   runtime_SGD=runtime_SGD, learn_rate=learn_rate))
    #     }
    #   } else {p = 0}
    # }
    
    
    # update pars_1 to pars_new
    pars_i = pars_new
    c = c + 1
    pars = rbind(pars, t(pars_i))
    times = c(times, Sys.time()-t_i)
    
    if (print_info==TRUE){
      print(c)
      cat("pars_new = ", pars_i, "\n")
      #cat("stats_diff = \n", round(diff,1), "\n")
      cat("moving average of stats_diff_abs = ", sum(tail(stats_diff_abs,20)), "\n", sep="")
    }
  }
  
  runtime_SGD = Sys.time() - t_i
  
  return(list(pars=pars, stats_diff=stats_diff, stats_diff_abs=stats_diff_abs, D=D_info, times=times))
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
              D=D, M=M, MD=MD, runtime_tot=runtime_tot, learn_rate=learn_rate))
}



find_ic = function(N_max, stats_obs, n_gen, n_treesXgen, n_stats, n_pars, pars_intervals, max_attempts){
  # matrices whcih will contain respective pars and mean(stats)
  stats = matrix(nrow=n_gen, ncol=n_stats)
  pars = matrix(nrow=n_gen, ncol=n_pars)
  
  for (i in 1:n_gen){
    
    attempts_exceeded = FALSE
    
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
      res = create_tree_mat_fixed_N(N_max, pars[i,1], pars[i,2:ncol(pars)], attempt, max_attempts, max_t)
      tree_mat = res$tree_mat
      att = res$attempt
      
      if (att >= max_attempts) {
        attempts_exceeded = TRUE
        break  # Exit the inner loop
      }
      
      # get extant tree
      tree_extant = DDD::L2phylo(unname(tree_mat), dropextinct=TRUE)
      # computing primary and LTT stats on observed tree
      stats_single = stats_PDM(tree_extant)
      # append to stats
      stats_j = rbind(stats_j, stats_single)
    }
    
    # append mean to general stats
    stats[i,] = colMedians(stats_j, na.rm=TRUE)
    
    cat("find_ic(): gen ", i, " done\n", sep="")
  }
  
  # find abs() of difference in stats
  stat_diff = abs(sweep(stats, 1, stats_obs))
  # calc row sums
  sums = rowSums(stat_diff)
  # find minimum
  id_min = which(sums == min(sums))
  # set ic to pars[id_max,]
  ic = pars[id_min,]
  
  return(ic)
}



#################################################################################################

create_tree_mat_fixed_N <- function(N_max, alphas, betas, attempt, max_attempts, max_t, print_info=FALSE) {
  init = init_tree_mat()
  tree_mat = init$tree_mat
  tips_at_t = init$tips_at_t
  
  counter = init$counter
  # set t = 0, present time tp, total branch length
  t = 0
  left = c(1)
  right = c(2)
  
  # Set counter for number of species
  n_count = 2
  
  # Gillespie algorithm to simulate birth-death process
  while (n_count <= N_max) {
    
    if (n_count>2){
      tree_mat_virt = tree_mat
      tree_mat_virt[,"birth.time"] = t - tree_mat_virt[,"birth.time"]
      #tree_mat_virt[,"death.time"] = t - tree_mat_virt[,"death.time"]
      tree_mat_virt[is.na(tree_mat_virt[,"death.time"]),"death.time"] = -1
      #tree_mat_virt = tree_mat_virt[1:(nrow(tree_mat_virt)-1),]
      # need to make calc_PDM faster
      PDM_full = calc_PDM(L2phylo(unname(tree_mat_virt), dropextinct=TRUE))
    } else
    {PDM_full = matrix(data=0,nrow=2,ncol=2)
    colnames(PDM_full)=c("t1","t2")}
    
    pdm = calc_avg_phylodiv(PDM_full)
    avg_phylodists = pdm$avg_phylodists
    tip_ids = pdm$tip_ids
    
    #tip_map = match(tip_ids, tips_at_t)
    
    tryCatch(
      {t_sample = rexp.var(1, function(tm) lambda_tot(tm, n_count, betas, avg_phylodists)+
                          mu_tot(tm, n_count, alphas, avg_phylodists), 
                          now=0, tMax=Inf, shape=NULL, TS=0)},
      error=function(e){
        return(c())
      }
    )
    
    
    # set new time
    t = t + t_sample
    
    if (t > max_t){
      attempt=max_attempts
      return(list(tree_mat=tree_mat, attempt=attempt, t_tot=t))
    }
    #print(t)
    #cat("N = ", n_count, "\n")
    #cat("mu_tot = ", mu_tot(t,n_count,alphas,avg_phylodists), "\n")
    #cat("lam_tot = ", lambda_tot(t,n_count,betas,avg_phylodists), "\n")
    
    # sample event type
    sampling_probs = 
      #constant part
      rep(exp(alphas[1]+alphas[length(alphas)]*n_count) + 
        exp(betas[1]+betas[length(betas)]*n_count), n_count) + 
      # lineage-varying part
        exp(alphas[2]*avg_phylodists) + exp(betas[2]*avg_phylodists)
    
    sampling_probs = sampling_probs/sum(sampling_probs)
    
    
    # allocate event to lineage by sampling from tips_at_t
    if (length(tips_at_t) != 0 & length(left) != 0 & length(right) != 0) {
      
      tryCatch(
        {node = sample(tip_ids, size=1, prob=unname(sampling_probs))},
        error = function(e){return(c())}
      )
      
    } else { 
      if (print_info==TRUE){
        cat("Generation interrupted as lineages have become extinct at time ", t, "\n")
        cat("Attempt n°", attempt, "\n")
      }
      if (attempt >= max_attempts) {
        cat("Max number of attempts reached: ", max_attempts, "\n")
        return(list(tree_mat=tree_mat, attempt=attempt, t_tot=t))
      } else {
        sub_result <- create_tree_mat_fixed_N(N_max, alphas, betas, attempt+1, max_attempts, max_t)
        sub_tree_mat <- sub_result$tree_mat
        sub_attempt <- sub_result$attempt
        
        if (!is.null(sub_tree_mat)) {
          tree_mat <- sub_tree_mat
        }
        attempt <- sub_attempt
        
        return(list(tree_mat=tree_mat, attempt=attempt, t_tot=t))
      }
    }
    
    # find the extinction rate percentage mu_tot / (lam_tot + mu_tot)
    ext_perc = mu_tot(0, n_count, alphas, avg_phylodists) /
      (lambda_tot(0, n_count, betas, avg_phylodists) + mu_tot(0, n_count, alphas, avg_phylodists))
    p = runif(1, min=0, max=1)
    
    # if event speciation
    if (p > ext_perc){
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
    
    n_count = length(tips_at_t)
    
    #cat("N = ", n_count, "\n")
    #cat("mu_tot = ", n_count*mu(mu_0), "\n")
    #cat("lambda_tot = ", lambda_tot(0, n_count, betas, avg_phylodists), "\n")
    #cat("mu_% = ", ext_perc, "\n")
    
  }
  
  
  # put matrix in correct form by setting tp = 0 and ti = tp
  tree_mat[,"birth.time"] = t - tree_mat[,"birth.time"]
  tree_mat[,"death.time"] = t - tree_mat[,"death.time"]
  
  # set NAs = -1 for extant species
  tree_mat[is.na(tree_mat[,"death.time"]),"death.time"] = -1
  
  tree_mat = tree_mat[1:(nrow(tree_mat)-1),]
  
  return(list(tree_mat=tree_mat, attempt=attempt, t_tot=t))
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
lambda_i <- function(id, t, N, betas, PDM) {
  sp = exp(betas[1] + betas[length(betas)] * N + betas[2]*calc_avg_phylodiv(id, PDM))
  return(sp)
}

lambda_tot <- function(t, N, betas, avg_phylodists){
  # Isolating the N  b[1]*b_N*N terms and the pdm terms
  l_tot = N * exp(betas[1]+betas[length(betas)]*N) * sum(exp(betas[2]*(avg_phylodists+(2*t))))
  
  return(l_tot)
}


mu_tot <- function(t, N, alphas, avg_phylodists){
  # Isolating the N  b[1]*b_N*N terms and the pdm terms
  m_tot = N * exp(alphas[1]+alphas[length(alphas)]*N) * sum(exp(alphas[2]*(avg_phylodists+(2*t))))
  
  return(m_tot)
}


calc_avg_phylodiv = function(PDM_full){
  # find tip ids
  tip_ids = colnames(PDM_full)
  tip_ids = as.numeric(gsub("t", "", tip_ids))
  # find means of columns excluding diagonal value (which is 0)
  avg_phylodists = colSums(PDM_full)/(ncol(PDM_full)-1)
  
  return(list(tip_ids=tip_ids, avg_phylodists=avg_phylodists))
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


calc_D_stats = function(n_stats, n_pars, pars_i, n_trees_D, max_attempts, max_t, print_info=TRUE){
  
  #pars_i = true_pars
  eps = pars_i/10
  
  s = matrix(nrow=0,ncol=n_stats)
  # generate original trees
  for (i in 1:n_trees_D){
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
    if (print_info==TRUE){
      print(i)
    }
  }
  stats_orig = colMedians(s)
  
  stats_grad = matrix(nrow=n_pars, ncol=n_stats)
  # generate trees with tweaked pars for each par
  for (p in 1:n_pars){
    #print(p)
    pars_eps = pars_i
    pars_eps[p] = pars_eps[p] + eps[p]
    
    s1 = matrix(nrow=0, ncol=n_stats)
    
    for (i in 1:n_trees_D){
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
      
      s1 = rbind(s1, stats_obs_pdm)
      if (print_info==TRUE){
        cat(p," ", i, "\n")
      }
    }
    # take average
    stats_avg = colMedians(s1)
    # populate row of matrix
    
    stats_grad[p,] = (stats_avg - stats_orig)/eps[p]
  }
  
  return(stats_grad)
}



calc_PDM <- function(tree) {
  terminals <- tree$tip.label
  n <- length(terminals)
  dist_matrix <- matrix(0, nrow = n, ncol = n)
  
  phylo_dist <- cophenetic(tree)  # Compute cophenetic distances once
  
  for (i in 1:(n - 1)) {
    dist_matrix[i, (i + 1):n] <- phylo_dist[terminals[i], terminals[(i + 1):n]]
  }
  
  dist_matrix <- t(dist_matrix) + dist_matrix  # Make it symmetric
  
  df <- as.data.frame(dist_matrix)
  colnames(df) <- terminals
  rownames(df) <- terminals
  
  return(df)
}

reorder_PDM <- function(matrix) {
  n <- ncol(matrix)
  
  # Find the smallest positive values in each column
  smallest_vals <- apply(matrix, 2, function(col) min(col[col > 0]))
  
  # Create an order based on the smallest values
  col_order <- order(smallest_vals)
  
  # Reorder the matrix's columns and rows according to this order
  reordered_matrix <- matrix[, col_order]
  reordered_matrix <- reordered_matrix[col_order, ]
  
  # set the bottom half to 0
  reordered_matrix[lower.tri(reordered_matrix)] <- 0
  
  return(reordered_matrix)
}


vectorize_PDM <- function(matrix) {
  n <- nrow(matrix)
  
  # Initialize the vector
  vector_length <- n * (n - 1) / 2
  vectorized_values <- numeric(vector_length)
  index <- 1
  
  for (i in 1:(n - 1)) {
    for (j in (i + 1):n) {
      vectorized_values[index] <- matrix[i, j]
      index <- index + 1
    }
  }
  return(vectorized_values)
}


stats_PDM <- function(tree) {
  PDM <- calc_PDM(tree)
  
  # Get the lower triangular part of the distance matrix (excluding diagonal)
  lower_triangular <- PDM[lower.tri(PDM, diag = FALSE)]
  
  # Sort the distances in ascending order
  sorted_distances <- sort(lower_triangular)
  
  return(sorted_distances)
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