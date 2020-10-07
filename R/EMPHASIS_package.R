### EMPHASIS functions
emphasisLD <- function(phylo,
                     init_par,
                     model,
                     em_tol=0.25,
                     sample_size_tol=0.005,
                     burnin_sample_size=200,
                     pilot_sample_size=c(200,600),
                     burnin_iterations = 20,
                     parallel=TRUE){

  brts = ape::branching.times(phylo)
  msg1 = paste("Initializing emphasis...")
  msg2 = paste("Age of the tree: ",max(brts))
  msg3 = paste("Number of speciations: ",length(brts))
  msg4 = paste("Diversification model to fit:",model)
  msg5 = "######################################"
  cat(msg1,msg2,msg3,msg4,msg5,sep="\n")
  
  cat( "Performing Phase 1: burn-in",sep= "\n")
  mc = mcEM(phylo, 
            pars = init_par, 
            sample_size = burnin_sample_size,
            model = model,
            parallel = parallel,
            print_process = FALSE,
            tol = em_tol,
            burnin = burnin_iterations)
  
  M = mc$mcem
  pars = c(mean(tail(M$par1,n = nrow(M)/2)),
           mean(tail(M$par2,n = nrow(M)/2)),
           mean(tail(M$par3,n = nrow(M)/2)),
           mean(tail(M$par4,n = nrow(M)/2)))
  
  cat("\n",msg5,sep="\n")
  cat( "Phase 2: Assesing required MC sampling size \n")

  for(i in 1:length(pilot_sample_size)){
    cat(paste("\n Sampling size: ",as.character(pilot_sample_size[i]),"\n"))
    mc = mcEM(brts = brts,
              pars = pars,
              sample_size = pilot_sample_size[i],
              model = model,
              soc = soc,
              parallel = parallel,
              print_process = FALSE,
              tol = em_tol,
              burnin = 10)
    ta = tail(mc$mcem,n = nrow(M)/2)
    pars = c(mean(ta$par1),mean(ta$par2),mean(ta$par3),mean(ta$par4))
    M = rbind(M,mc$mcem)
  }
  n.r = get_required_sampling_size(M[-(1:burnin_iterations),],tol = sample_size_tol)
  sample_size = max(pilot_sample_size+2,n.r)
  n.r_old = -1
  j = 1
  while(n.r_old < n.r){
    msg6 = paste0("Required sampling size: ",n.r)
    msg7 = paste0("Phase 3: Performing metaiteration: ",j)
    cat("\n",msg5,msg7,msg6,sep="\n")
    mc = mcEM(brts = brts,
              pars = pars,
              sample_size = sample_size,
              model = model,
              soc = soc,
              parallel = parallel,
              print_process = FALSE,
              tol = em_tol,
              burnin = 2)
    M <- rbind(M,mc$mcem)
    n.r_old = n.r
    j = j+1
    n.r = get_required_sampling_size(M[-(1:burnin_iterations),],tol = sample_size_tol)
    pars = as.numeric(colMeans(mc$mcem)[1:4])
    sample_size = n.r
  }

  cat(pars)
  return(list(pars=pars,MCEM=M))
}

mcEM <- function(phylo, 
                 diversification_model, 
                 sample_size, 
                 tol=0.01,
                 burnin=20,
                 vebose=TRUE,
                 parallel=TRUE,
                 cores=(parallel::detectCores()-1)){
  mcem = NULL
  sde = 10; i=0
  times = NULL
  PARS = NULL
  while(sde > tol){
    i = i+1
    results = em_r(phylo,
                   diversification_model,
                   sample_size = sample_size, 
                   no_cores = cores, 
                   parallel = parallel)
    
    pars = results$estimates
    PARS = rbind(PARS,pars)
    log_fhat = c(log_fhat, results$log_fhat)

    if (i > burnin + 2) {
      
      tail_mcem = tail(log_fhat,length(log_fhat)-burnin)
      sde = sd(tail_mcem) / sqrt(length(tail_mcem))
      if (verbose) {
        msg <- paste("Iteration:", i, " SE of the loglikelihood: ", sde)
        cat("\r", msg)
      }
    } else {
      if (verbose) {
         msg <- paste("Performing burn-in, should take short time ")
         cat("\r", msg)
      }
    }
  }
  rownames(PARS) = NULL
  colnames(PARS) = paste0("par",1:ncol(PARS))
  
  return(list(PARS = PARS,log_fhat=log_fhat))
}

##############################

em_r <- function(phylo,
                 diversification_model,
                 sample_size = sample_size, 
                 no_cores = cores, 
                 parallel = parallel){
  ST = MC_augmentation(phylo,
                diversification_model,
                sample_size = sample_size, 
                no_cores = cores, 
                parallel = parallel)
  w = get_weights(ST = ST,
                  diversification_model = diversification_model)
  if(max(w$weight)==0){
    print(st)
    stop("Only zero likelihood trees, maybe there is underflow")
  }
  st = list(trees=ST$trees,weights=w$weights)
  M = M_step(st = st, init_par = pars, model = diversification_model$model)
  return(list(estimates=M$po$value,log_fhat=w$log_fhat))
}

####### E-step 
MC_augmentation <- function(phylo,
                     diversification_model,
                     sample_size,
                     no_cores=2,
                     seed=0,
                     parallel=TRUE){
  if(seed>0) set.seed(seed)

  if(!parallel){
    st =  lapply(1:sample_size,function(i){sample_tree(diversification_model = diversification_model,phylo = phylo)} )
  }else{
    st = mclapply(1:sample_size,function(i){sample_tree(diversification_model = diversification_model,phylo = phylo)},mc.cores = no_cores)
  }
  
  return(st)

}

get_weights <- function(ST,diversification_model){
  pars = diversification_model$pars
  model = diversification_model$model
  logf = sapply(ST,loglik.tree(model), pars=pars)
  logg = sapply(ST,sampling_prob, pars=pars,model=model)
  log_weights = logf-logg
  prop_const = max(log_weights)
  log_weights_norm = log_weights - prop_const
  w = exp(log_weights_norm)
  log_fhat = log(mean(w)) + prop_const
  return(list(weights=w,log_fhat=log_fhat))
}


##############################
####### M-step 

M_step <-function(st,init_par,model,reltol=0.001){
  
  time0 = proc.time()
  
  sub_st = get_contributing_trees(st)
  
  po = subplex(par = init_par, fn = Q_approx,st = sub_st, loglik = get(paste0("loglik.tree.", model)), hessian = FALSE,control=list(reltol=reltol))

  M_time = get.time(time0)
  return(list(po=po,M_time=M_time))
}


get_contributing_trees <- function(st, exclude_proportion_trees=0){
  w = st$weights/sum(st$weights)
  contributing_trees = (w > exclude_proportion_trees)
  sub_trees = st$trees[contributing_trees]
  effective_sample_size = sum(contributing_trees)
  sub_st = list(trees = sub_trees, weights = st$weights[contributing_trees])
  return(sub_st)
}

Q_approx = function(pars,st,loglik){
  
  l = sapply(st$trees, loglik, pars=pars)
  w = st$weights
  Q = -sum(l*w)
  return(Q)
}

####
