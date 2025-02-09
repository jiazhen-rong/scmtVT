#' This function is for Coverage-specific False Positives Estimation,
#' calculating the number of expected false positives in the significant cells.
#'
#' @param bin_size The bin size of coverage. Default is 20 counts per bin.
#' @param X_sub alternate allele count matrix.
#' @param N_sub total variant count matrix.
#' @param alpha parameter of beta prior, estimated from EM.
#' @param beta  parameter of beta prior, estimated from EM.
#' @param delta_ij_hat  parameter representing non-zero distribution component, estimated from EM.
#' @param FDR_alpha A numeric value representing FDR threshold.
#' @param gc  name of background cluster, used for saving the file.
#' @param save_path path to save the diagnostic plots. Default setting will store in "zibb_fit" folder under current working directory.
#' @param cell_label A vector containing cell identities based on scRNA-Seq and standard Seurat annotation
#' @param seu A seurat object generated from standard scRNA-Seq processing.
#' @param voi A list of Variant of Interest after basic filtering.
#' @param output TRUE or FALSE, representing whether to write significant cells of each cluster.
#' @param plot_umap TRUE or FALSE, representing whether to plot significant cells of each cluster on UMAP.
#' @param plot_diagnostic TRUE or FALSE, representing whether to plot diagnostic plot.
#' @param verbose TRUE or FALSE, representing whether to print out detailed log information.
#'
#' @return A list containing:
#' \itemize{
#'   \item sig_cell_df Dataframe containing number of significant cells in each variant.
#'   \item FP_df Dataframe containing number of expected false positives in each variant.
#' }
#' @export
#'
#' @examples
#' pval_df=ZIBB_test(X_sub,N_sub,alpha,beta,delta_ij_hat,gc="Normal",cell_label=cell_label,seu=seu,
#'           output=T,plot_umap=T,save_path="example/output/zibb_fit/plots/")
#'
CoverageFP <-function(X_sub,N_sub,alpha,beta,delta_ij_hat,bin_size = 20,voi,gc,FDR_alpha=0.1,cell_label,output=T,save_path=NULL){
  if(is.null(save_path)){
    save_path="zibb_fit/test_result/false_positive_rates/"
  }else{
    dir.create(save_path,recursive = T)
  }
  # variant-specific contamination proportion from background/normal cells for variants of interest
  pjN = (alpha + colSums(delta_ij_hat*X_sub,na.rm=TRUE))/(alpha + beta + colSums(delta_ij_hat*N_sub,na.rm=TRUE))
  # Coverage-Specific FDR
  sig_cell_df = c()
  FP_df = c()
  for(id in 1:length(voi)){
    message(paste0("Estimating False Positive number for variant ",id,": " ,voi[id],"..."))
    x=X[voi[id],]
    n=N[voi[id],]
    af=x/(n+0.1)
    # normal
    x0 = x[cell_label==gc]
    n0 = n[cell_label==gc]

    ## record rejected cells in normal
    # Binomial Test
    pnull=pjN[id]
    u1 = 1-pbinom(x0-1, n0, pnull)
    # FDR selection process in background/normal cells
    u1_order=order(u1)
    #below_thresh = u1[u1_order]<alpha*c(1:length(u1))/length(u1)
    below_thresh = u1[u1_order]<FDR_alpha*c(1:length(u1))/length(cell_label != gc)
    reject = sum(below_thresh)
    rejected_0=u1_order[which(below_thresh)] # rejected cells index in n0

    # test within each cluster
    sig_cell_num = c()
    total_FP=c()
    for (cluster in levels(cell_label)){
      if(cluster != "test"){ # skip normal cluster
        sel1 = which(cell_label==cluster)
        x1 = x[cell_label==cluster]
        n1= n[cell_label==cluster]
        # Binomial Test
        pnull=pjN[id]
        u1 = 1-pbinom(x1-1, n1, pnull)
        # FDR selection process
        u1_order=order(u1)
        below_thresh = u1[u1_order]<FDR_alpha*c(1:length(u1))/length(cell_label != "Normal")
        reject = sum(below_thresh)
        rejected=u1_order[which(below_thresh)]
        rejected_ind = sel1[rejected]
        # Add in UMAP information
        sig_cells = rep(0,dim(N)[2])
        sig_cells[rejected_ind] = 1 # 1 for significant, 0 for non-significant
        sig_cell_num <- append(sig_cell_num,length(rejected_ind))

        # calculate bin-wise rejections and rejections
        total_bin_num = ceiling(max(N[voi[id],])/bin_size)
        fpr = rep(NA,total_bin_num) # false positive rate
        fp_n = rep(NA,total_bin_num)
        fdr = rep(NA,total_bin_num) # false discovery rate
        FP_n = rep(NA,total_bin_num)
        for(bin_n in 1:total_bin_num){
          # fpr in empirical null distribution
          fpr[bin_n] = sum((n0[rejected_0] > bin_size *(bin_n-1)) &
                           (n0[rejected_0] <= bin_size *(bin_n))) /
                           (sum((n0 > bin_size *(bin_n-1)) & (n0 <= bin_size * bin_n)) + 0.001)
          fp_n[bin_n] = fpr[bin_n] *sum((n1 > bin_size *(bin_n-1)) & (n1 <= bin_size * bin_n))
          fdr[bin_n] = min(1,fp_n[bin_n]/(sum((n1[rejected] > bin_size *(bin_n-1)) &
                                              (n1[rejected] <= bin_size *(bin_n)))+0.001))
          FP_n[bin_n] = fdr[bin_n] * (sum((n1[rejected] > bin_size *(bin_n-1)) &
                                          (n1[rejected] <= bin_size *(bin_n)))+0.001)
        }
        total_FP = append(total_FP,sum(FP_n))
      }
    }
    FP_df <- rbind(FP_df,total_FP)
    sig_cell_df <- rbind(sig_cell_df,sig_cell_num)
  }
  sig_cell_df <- as.data.frame(sig_cell_df)
  colnames(sig_cell_df) <- levels(cell_label)
  rownames(sig_cell_df) <- voi

  FP_df <- as.data.frame(FP_df)
  colnames(FP_df) <- levels(cell_label)
  rownames(FP_df) <- voi
  if(output==T){
    write.csv(sig_cell_df,
          file=paste0(save_path,"/sig_cell_number.txt"))
    write.csv(FP_df,
          file=paste0(save_path,"/expected_false_positive_number.txt"))
  }
  return(list(sig_cell_df=sig_cell_df,FP_df=FP_df))
}
