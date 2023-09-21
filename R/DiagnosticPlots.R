#' Create diagnostic plot of the ZIBB EM parameter inference process
#'
#' @param X_sub alternate allele count matrix.
#' @param N_sub total variant count matrix.
#' @param alpha initial value of alpha
#' @param beta  initial value of beta
#' @param delta_gc initial value of delta_gc
#' @param delta_ij_hat initial value of delta_ij_hat, set to be the same as delta_gc
#' @param gc  name of background cluster, used for saving the file.
#' @param save_path path to save the diagnostic plots. Default setting will store in "zibb_fit" folder under current working directory.
#' @param delta_gc_list list of parameters representing non-zero distribution component in each iteration,  estimated from EM process.
#' @param alpha_list list of parameters of beta prior in each iteration, estimated from EM process all iterations.
#' @param beta_list list of parameters of beta prior in each iteration, estimated from EM process all iterations.
#' @param LL_list list of log likelihood in each iteration.
#' @param delta_ij_hat matrix representing missing variable Zij
#' @param width width of saved pdf, unit: inches
#' @param height height of saved pdf, unit: inches
#'
#' @return A diagnostic PDF.
#' @export
#'
#' @examples
#' ZIBB_fit_diagnostic_plot(X_sub,N_sub,delta_gc=0.5,alpha=1,beta=1,delta_ij_hat=0.5,gc="normal",
#'                          delta_gc_list,LL_list,alpha_list,beta_list,save_path="example/output/zibb_fit/",
#'                          width=10,height=10)
#'
ZIBB_fit_diagnostic_plot <- function(X_sub,N_sub,delta_gc=0.5,alpha=1,beta=1,delta_ij_hat=0.5,gc="normal",
                                     delta_gc_list,LL_list,alpha_list,beta_list,save_path=NULL,width=10,height=10){
  # plot diagnostic plot & save them
  i=length(alpha_list)
  init_LL = loglik(X_sub,N_sub,delta_gc=delta_gc,alpha=alpha,beta=beta,delta_ij_hat=delta_ij_hat)
  if(is.null(save_path)){
    pdf(paste0("zibb_fit/start_a_",1,"_b_",1,"_del_",0.5,"_",gc,".pdf"),width=width,height=height)
  }else{
    dir.create(save_path,recursive = T)
    pdf(paste0(save_path,"/start_a_",alpha,"_b_",beta,"_del_",delta_gc,"_",gc,".pdf"),width=width,height=height)
  }
  par(mfcol= c( 2, 2 ))
  plot(1:(i+1),c(delta_gc,delta_gc_list),type="l",main="Change of delta_gc",xlab="Iterations",ylab="delta_gc")
  plot(1:(i+1),c(init_LL,LL_list),type="l",main="Change of Log Likelihood",xlab="Iterations",ylab="loglik")
  plot(1:(i+1),c(alpha,alpha_list),type="l",main="Change of alpha",xlab="Iterations",ylab="alpha")
  plot(1:(i+1),c(beta,beta_list),type="l",main="Change of beta",xlab="Iterations",ylab="beta")

  par(mfcol= c( 2, 2 ))
  plot(1:i,delta_gc_list,type="l",main="Change of delta_gc (Exclude start)",xlab="Iterations",ylab="delta_gc")
  plot(1:i,LL_list,type="l",main="Change of Log Likelihood (Exclude start)",xlab="Iterations",ylab="loglik")
  plot(1:i,alpha_list,type="l",main="Change of alpha (Exclude start)",xlab="Iterations",ylab="alpha")
  plot(1:i,beta_list,type="l",main="Change of beta (Exclude start)",xlab="Iterations",ylab="beta")

  dev.off()
  message("Diagnostic Plot of ZIBB Process created!")
}
