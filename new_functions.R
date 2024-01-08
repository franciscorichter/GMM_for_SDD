library(parallel)

# GRADIENT DESCENT FOR EVOLUTIONARY MODEL WITH CONSTANT D MATRIX
gradDescentEvolutionModel <- function(N_max, stats_obs, n_stats, n_pars, learn_rate, patience,
                                      iters, initial_params, n_trees_D, n_trees_sgd, max_t, max_attempts,print_info=TRUE) {
  iteration_counter <- 0
  patience_counter <- 0
  pars_matrix <- matrix(nrow=0, ncol=length(initial_params))
  pars_matrix <- rbind(pars_matrix, t(initial_params))
  updated_params <- initial_params
  stats_diff_matrix <- matrix(nrow=0, ncol=n_stats)
  stats_diff_abs_vector <- c()
  execution_times <- c()
  
  previous_params <- updated_params
  
  # Initial scaling
  scaling <- abs(updated_params)
  D_matrix_info <- list()
  stats_diff_history <- NULL
  
  convergence_threshold <- 1e-3  
  
  while (iteration_counter < iters) {
    iteration_start_time <- Sys.time()
    
    # Update D matrix periodically
  
    if (iteration_counter == 0 | isStatsDiffIncreasing(stats_diff_history)) {
      if (print_info) cat("Calculating D matrix...")
      D_matrix <- calculateDMatrix(n_stats = n_stats, 
                                   n_pars = n_pars, 
                                   initial_params = initial_params,
                                   n_trees_D =  n_trees_D,
                                   max_attempts =  max_attempts,
                                   max_t =  max_t,
                                   use_parallel = FALSE)
      if (length(D_matrix) == 0) { D_matrix <- D_matrix_info[[length(D_matrix_info)]] }
      
      MD_matrix <- D_matrix / rowSums(abs(D_matrix))
      D_matrix_info <- c(D_matrix_info, list(D_matrix=D_matrix,iteration=iteration_counter))
    }
    
    generated_stats <- generateStatistics(N_max, updated_params, n_trees_sgd, max_attempts, max_t)
    
    # Check for errors in tree generation
    if (is.null(generated_stats)) {
      return(list(pars=pars_matrix, stats_diff=stats_diff_matrix, stats_diff_abs=stats_diff_abs_vector, D=D_matrix_info, times=execution_times))
    }
    
    # Calculate difference in statistics
    stats_diff <- colMedians(generated_stats) - stats_obs
    
    # Update parameters
    gradient_adjustment <- learn_rate * (MD_matrix %*% (as.matrix(stats_diff))) * scaling
    updated_params <- initial_params - gradient_adjustment
    
    # Update condition
    parameter_change <- updated_params - initial_params
    stats_diff_matrix <- rbind(stats_diff_matrix, stats_diff)
    stats_diff_abs_vector <- c(stats_diff_abs_vector, sum(abs(stats_diff)))
    stats_diff_history <- c(stats_diff_history, sum(abs(stats_diff)))

    
    
    # Update initial parameters for next iteration
    pars_matrix <- rbind(pars_matrix, t(initial_params))
    initial_params <- updated_params
    iteration_counter <- iteration_counter + 1
    execution_times <- c(execution_times, Sys.time() - iteration_start_time)
    
    previous_params <- updated_params
    
    # Logging
    if (print_info) {
      print(iteration_counter)
      cat("Updated Params =", initial_params, "\n")
      cat("Moving average of stats_diff_abs =", mean(tail(stats_diff_abs_vector, 20)), "\n")
    }
    
    if (sum(abs(updated_params - previous_params)) < convergence_threshold) {
    #  break  # Exit the loop
    }
    
  }
  pars_matrix <- rbind(pars_matrix, t(initial_params))
  total_runtime <- Sys.time() - iteration_start_time
  return(list(pars=pars_matrix, stats_diff=stats_diff_matrix, stats_diff_abs=stats_diff_abs_vector, D=D_matrix_info, times=execution_times))
}

isStatsDiffIncreasing <- function(history) {
  if (length(history) < 10) return(FALSE)
  return(mean(tail(history, 10)) > mean(tail(history, 9)) & 
           mean(tail(history, 9)) > mean(tail(history, 8))& 
           mean(tail(history, 8)) > mean(tail(history, 7))) 
}


generateStatistics <- function(N_max, params, n_trees_sgd, max_attempts, max_t) {
  stats_matrix <- matrix(nrow=0, ncol=n_stats)
  
  for (i in 1:n_trees_sgd) {
    result <- create_tree_mat_fixed_N(N_max, params[1:3], params[4:length(params)], 0, max_attempts, max_t)
    if (length(result) == 0 || result$attempt >= max_attempts) {
      return(NULL)
    }
    
    tree_mat <- result$tree_mat
    tree <- L2phylo(unname(tree_mat), dropextinct=TRUE)
    stats_matrix <- rbind(stats_matrix, stats_PDM(tree))
  }
  
  return(stats_matrix)
}



calc_stats_for_perturbed_params <- function(pars, N_max, n_stats, n_pars, n_trees_D, max_attempts, max_t) {
  s = matrix(nrow=0, ncol=n_stats)
  for (i in 1:n_trees_D) {
    result = create_tree_mat_fixed_N(N_max, pars[1:3], pars[4:n_pars], 0, max_attempts, max_t)
    if (length(result) == 0 || result$attempt >= max_attempts) {
      return(NULL)
    }
    tree_mat = result$tree_mat
    tree_extant = L2phylo(unname(tree_mat), dropextinct=TRUE)
    s = rbind(s, stats_PDM(tree_extant))
  }
  return(colMedians(s))
}





calc_D_stats <- function(n_stats, n_pars, pars_i, n_trees_D, max_attempts, max_t, print_info = TRUE) {
  eps <- pars_i / 10
  
  # Helper function to generate stats
  generate_stats <- function(pars) {
    s <- matrix(nrow = 0, ncol = n_stats)
    for (i in 1:n_trees_D) {
      result <- create_tree_mat_fixed_N(N_max, pars[1:3], pars[4:n_pars], 0, max_attempts, max_t)
      if (is.null(result) || result$attempt >= max_attempts) {
        return(NULL)
      }
      tree_mat <- result$tree_mat
      tree_extant <- L2phylo(unname(tree_mat), dropextinct = TRUE)
      s <- rbind(s, stats_PDM(tree_extant))
    }
    colMedians(s)
  }
  
  # Calculate original stats
  stats_orig <- generate_stats(pars_i)
  if (is.null(stats_orig)) {
    stop("Failed to generate original stats.")
  }
  
  # Initialize gradient matrix
  stats_grad <- matrix(nrow = n_pars, ncol = n_stats)
  
  # Calculate gradients
  for (p in 1:n_pars) {
    pars_eps <- pars_i
    pars_eps[p] <- pars_eps[p] + eps[p]
    stats_eps <- generate_stats(pars_eps)
    if (is.null(stats_eps)) {
      stop(paste("Failed to generate stats for parameter", p))
    }
    stats_grad[p, ] <- (stats_eps - stats_orig) / eps[p]
    if (print_info) {
      cat("Parameter", p, "done\n")
    }
  }
  
  stats_grad
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


calc_avg_phylodiv = function(PDM_full){
  # find tip ids
  tip_ids = colnames(PDM_full)
  tip_ids = as.numeric(gsub("t", "", tip_ids))
  # find means of columns excluding diagonal value (which is 0)
  avg_phylodists = colSums(PDM_full)/(ncol(PDM_full)-1)
  
  return(list(tip_ids=tip_ids, avg_phylodists=avg_phylodists))
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

stats_PDM <- function(tree) {
  PDM <- calc_PDM(tree)
  
  # Get the lower triangular part of the distance matrix (excluding diagonal)
  lower_triangular <- PDM[lower.tri(PDM, diag = FALSE)]
  
  # Sort the distances in ascending order
  sorted_distances <- sort(lower_triangular)
  
  return(sorted_distances)
}

calculateDMatrix <- function(n_stats, n_pars, initial_params, n_trees_D, max_attempts, max_t, use_parallel = FALSE) {
 # if (use_parallel) {
    # Call the parallel version of the function
 #   stats_grad <- calc_D_stats_parallel(n_stats, n_pars, initial_params, n_trees_D, max_attempts, max_t, print_info=FALSE)
#  } else {
    # Call the non-parallel version of the function
    stats_grad <- calc_D_stats(n_stats, n_pars, initial_params, n_trees_D, max_attempts, max_t, print_info=FALSE)
 # }
  return(stats_grad)
}
