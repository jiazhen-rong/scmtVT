#' Expectation Step (E-step) to calculate the expected value of missing variable Zij.
#'
#' @param X_sub alternate allele count matrix.
#' @param N_sub total variant count matrix.
#' @param alpha parameter of beta prior, estimated from M-step.
#' @param beta  parameter of beta prior, estimated from M-step.
#' @param delta_gc  parameter representing non-zero distribution component, estimated from M-step.
#'
#' @return delta_ij_hat, a matrix representing the expected value of each cell and variant's missing identity Zij.
#' @export
#'
#' @examples
#' delta_ij_hat = Estep(X_sub,N_sub,alpha,beta,delta_gc)
Estep <- function(X_sub,N_sub,alpha,beta,delta_gc){
  # If Xij = 0, could come from either zero distribution or beta-bin distribution. Xij = 0 for ~95% times.
  # calculate delta_ij for each variant (# this is way faster implementation)
  delta_ij_hat = sapply(1:dim(X_sub)[2],function(j){
    P_betabin= dbetabinom.ab(X_sub[,j], size=N_sub[,j], shape1=alpha, shape2=beta)
    P_betabin*delta_gc/(P_betabin*delta_gc+1-delta_gc)
  })
  
  #  If Xij > 0 case, delta_ij = Zij = 1, variant only from Beta-Bin distirbution
  delta_ij_hat[which(X_sub!=0,arr.ind = TRUE)] = 1
  # ~ 2% of Nij (total count) are 0s, setting these probability to NA
  delta_ij_hat[which(N_sub==0,arr.ind = TRUE)] = NA
  return(delta_ij_hat)
}

#' Negative log likelihood related with alpha and beta
#'
#' @param parm a list containing the rest parameters in the ZIBB model.
#' @param data alist containing all data in the ZIBB model.
#'
#' @return NLL, a numeric value representing negative log likelihood related with alpha and beta
#' @import VGAM
#' @export 
#' 
#' @examples
#' NLL=negative_loglik_ab_related(parm,data)
negative_loglik_ab_related  <-function(parm,data){
  alpha<-parm[1]; beta <-parm[2]
  # negative log likelihood
  NLL = -sum( data$delta_ij_hat * dbetabinom.ab(data$X_sub,
                                                size=data$N_sub,
                                                shape1=alpha, shape2=beta,log =TRUE),na.rm=TRUE)
  return(NLL)
}

#' Total log likelihood of the whole model and all data
#'
#' @param X_sub alternate allele count matrix.
#' @param N_sub total variant count matrix.
#' @param alpha parameter of beta prior, estimated from M-step.
#' @param beta  parameter of beta prior, estimated from M-step.
#' @param delta_gc  parameter representing non-zero distribution component, estimated from M-step.
#' @param delta_ij_hat a matrix representing the expected value of each cell and variant's missing identity Zij.
#'
#' @return LL, a numeric value representing total log likelihood 
#' @import VGAM
#' @export
#'
#' @examples
#' new_LL = loglik(X_sub,N_sub,delta_gc,alpha,beta,delta_ij_hat)
loglik <- function(X_sub,N_sub,delta_gc,alpha,beta,delta_ij_hat){
  LL = sum( c(as.matrix(delta_ij_hat)) * dbetabinom.ab(c(as.matrix(X_sub)),
                                                       size=c(as.matrix(N_sub)),
                                                       shape1=alpha, shape2=beta,log =TRUE),na.rm=TRUE) +
    sum(delta_ij_hat*log(delta_gc) + (1-delta_ij_hat)*log(1-delta_gc),na.rm=TRUE)
  return(LL)
}

#' Maximization Step in the EM, estimate each parameter based on expectation of Zij from E-step.
#'
#' @param X_sub alternate allele count matrix.
#' @param N_sub total variant count matrix.
#' @param alpha_old parameter of beta prior, from last iteration.
#' @param beta_old  parameter of beta prior, from last iteration.
#' @param opt_method  parameter of optimization method for alpha and beta. optim: optimize with R's optim pacakge.
#'
#' @return a list of parameters estimated from the M-step:
#' \itemize{
#'   \item delta_gc_new - parameter representing non-zero distribution component, estimated from M-step in each iteration.
#'   \item alpha_new - parameter of beta prior, estimated from M-step in each iteration.
#'   \item beta_new - parameter of beta prior, estimated from M-step in each iteration.
#' }
#' @export
#'
#' @examples
#' params = Mstep(X_sub,N_sub,delta_ij_hat,alpha,beta,opt_method=c("optim"))
Mstep <- function(X_sub,N_sub,delta_ij_hat,alpha_old,beta_old,opt_method=c("optim","MoM")){
  delta_gc_new = sum(delta_ij_hat,na.rm=TRUE)/(dim(X_sub)[1]*dim(X_sub)[2] - sum(N_sub==0)) # adjusted for Nij = 0's NA delta_ij values
  # optim to find alpha and beta -- way too slow
  if(opt_method == "optim"){ # R's optimization function
    #start.time <- Sys.time()
    param_fit = optim(c(alpha_old,beta_old), negative_loglik_ab_related,
                      data=data.frame(delta_ij_hat=c(as.matrix(delta_ij_hat)),X_sub=c(as.matrix(X_sub)),N_sub=c(as.matrix(N_sub))),
                      method = "L-BFGS-B",lower=c(1e-5,1e-5), upper=c(Inf,Inf))
    #end.time <- Sys.time()
    #time.taken <- end.time - start.time
    alpha_new = param_fit$par[1]
    beta_new = param_fit$par[2]
  } else if("NR"){ # Newton-Raptson for estimating alpha and beta
    
  } # else if("MoM"){ # methods of moments for alpha and beta
  #}
  return(list(delta_gc_new,alpha_new,beta_new))
}


#' The integrative function for expectation maximization
#'
#' @param X_sub alternate allele count matrix.
#' @param N_sub total variant count matrix.
#' @param alpha parameter of beta prior, estimated from M-step.
#' @param beta  parameter of beta prior, estimated from M-step.
#' @param delta_gc  parameter representing non-zero distribution component, estimated from M-step.
#' @param iterations  number of EM iterations.
#' @param stop_diff  threshold of stopping the EM process if log likelihood differentce is smaller than the threshold.
#' @param gc  name of background cluster, used for saving the file.
#' @param verbose whether to print out the detailed time and debugging information in every 1/10 iterations.
#' @param save_path path to save the fitted parameters. Default setting will store in "zibb_fit" folder under current working directory.
#'
#' @return A list of estimated parameters:
#' \itemize{
#'   \item delta_gc_list - list of parameters representing non-zero distribution component in each iteration,  estimated from EM process.
#'   \item alpha_list - list of parameters of beta prior in each iteration, estimated from EM process all iterations.
#'   \item beta_list - list of parameters of beta prior in each iteration, estimated from EM process all iterations.
#'   \item LL_list - list of log likelihood in each iteration.
#'   \item delta_ij_hat - matrix representing missing variable Zij 
#' }
#' 
#' @import VGAM
#' @export 
#' @examples
#' param_fit = EM_ZIBB(X_sub,N_sub,alpha=1,beta=1,delta_gc=0.5,iterations=100,stop_diff=1e-05,save_path=NULL,gc="normal",verbose=T)
EM_ZIBB <- function(X_sub,N_sub,alpha=1,beta=1,delta_gc=0.5,iterations=10,stop_diff=1e-05,gc="Normal",save_path=NULL,verbose=T){
  init_alpha=alpha;init_beta=beta;init_delta_gc=delta_gc;
  delta_gc_list = c()
  alpha_list = c()
  beta_list = c()
  LL_list = c()
  prev_param = list(0,0,0)
  params = list(alpha,beta,delta_gc)
  start.time0 <- Sys.time()
  for(i in 1:iterations){
    #print(i)
    diff = sum(unlist(prev_param) - unlist(params))
    if(abs(diff) < stop_diff){
      print("equal")
      break
    }
    #start.time <- Sys.time()
    # Estep to get 
    prev_param = params
    delta_ij_hat = Estep(X_sub,N_sub,alpha,beta,delta_gc)
    # Maximization Step in the EM, estimate each parameter based on expectation of Zij
    params = Mstep(X_sub,N_sub,delta_ij_hat,alpha,beta,opt_method=c("optim"))
    # save all parameters for diagnostic
    delta_gc=params[[1]]
    delta_gc_list = append(delta_gc_list,delta_gc)
    alpha=params[[2]]
    alpha_list = append(alpha_list,alpha)
    beta=params[[3]]
    beta_list = append(beta_list,beta)
    new_LL = loglik(X_sub,N_sub,delta_gc,alpha,beta,delta_ij_hat)
    LL_list = append( LL_list,new_LL)
    if((verbose==T) & (i%%round(iterations/10)==0)){
      #end.time <- Sys.time()
      #time.taken <- end.time - start.time
      end.time0 <- Sys.time()
      time.taken0 <- difftime(end.time0, start.time0, units='mins') #end.time0 - start.time0
      message(paste0("Iteration till ", i, "th takes ", format(time.taken0,digits=2)," s"))
    }
  }
  end.time0 <- Sys.time()
  time.taken0 <- difftime(end.time0, start.time0, units='mins') #end.time0 - start.time0
  message(paste0("Total EM time for ", i, " iterations take ",format(time.taken0,digits=2)," mins"))

  # save the estimated parameter for cluster c (gc)
  if(is.null(save_path)){
    dir.create("zibb_fit/")
    save(delta_gc_list,alpha_list,beta_list, LL_list,delta_ij_hat,
     file = paste0("zibb_fit/start_a_",init_alpha,"_b_",init_beta,"_del_",init_delta_gc,"_",gc,".RData"))
  }else{
    dir.create(save_path,recursive = T)
    save(delta_gc_list,alpha_list,beta_list, LL_list,delta_ij_hat,
         file = paste0(save_path,"/start_a_",init_alpha,"_b_",init_beta,"_del_",init_delta_gc,"_",gc,".RData"))
  }
  return(list(delta_gc_list=delta_gc_list,alpha_list=alpha_list,beta_list=beta_list,LL_list=LL_list,delta_ij_hat=delta_ij_hat))
}
