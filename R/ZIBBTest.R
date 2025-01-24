#' This function contains corresponding ZIBB test for variant of interest
#' and creates corresponding diagnostic plots.
#' The significant cells for each variant of interest will be labeled in the diagnostic plot and the UMAP.
#'
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
#' @return A matrix containing ZIBB test's p-values of each cell i for each given variant of interest j (cell x variant)
#' @import ggplot2 grid gridExtra RColorBrewer
#'
#' @export
#'
#' @examples
#' pval_df=ZIBB_test(X_sub,N_sub,alpha,beta,delta_ij_hat,gc="Normal",cell_label=cell_label,seu=seu,
#'           output=T,plot_umap=T,save_path="example/output/zibb_fit/plots/")
#'
ZIBB_test <- function(X_sub,N_sub,alpha,beta,delta_ij_hat,FDR_alpha=0.1,gc="Normal",cell_label=NULL,seu=NULL,voi=NULL,
                      output=T,plot_umap=T,save_path=NULL,plot_diagnostic=T,verbose=F,
                      umap_coord=NULL,celltype_color=NULL,celltype_cex=NULL,celltype_pch=NULL){
  if(is.null(save_path)){
    save_path="zibb_fit/plots/"
  }else{
    dir.create(save_path,recursive = T)
  }
  if(is.null(umap_coord)){
    if(is.null(seu@reductions$umap)){
      stop("No UMAP embeddings, please either running standard UMAP or input UMAP coordinates with umap_coord parameter.")
    }
    umap_coord=seu@reductions$umap@cell.embeddings
  }
  # checking for inputs
  if (!is.matrix(X_sub) || !is.matrix(N_sub)) stop("X_sub and N_sub must be matrices")
  if (dim(X_sub) != dim(N_sub)) stop("X_sub and N_sub must have the same dimensions")
  if (!is.factor(cell_label)) stop("cell_label must be a factor")
  if (!all(levels(cell_label) %in% unique(cell_label))) stop("Mismatch in cell_label levels")
  if (!is.null(seu) && !"Seurat" %in% class(seu)) stop("seu must be a Seurat object")

  # variant-specific contamination proportion from background/normal cells for variants of interest
  pjN = (alpha + colSums(delta_ij_hat*X_sub,na.rm=TRUE))/(alpha + beta + colSums(delta_ij_hat*N_sub,na.rm=TRUE))
  # dataframe of saving p-value of each given variant of interest
  pval_df = c()
  # Plot diagnostic plot for each variant
  for(id in 1:length(voi)){
    message(paste0("Testing for variant ",id,": " ,voi[id],"..."))
    diagnostic_plot_path = paste0(save_path,"/",voi[id], "_all_diagnostic_plot.pdf")
    umap_path = paste0(save_path,"/",voi[id], "_significant_cells_UMAP.pdf")

    x=X[voi[id],]
    n=N[voi[id],]
    af=x/(n+0.1) # add in a small value to avoid division by 0
    # normal
    x0 = x[cell_label==gc]
    n0 = n[cell_label==gc]

    # UMAP information

    plot_df <- as.data.frame(cbind(umap_coord, #seu@reductions$umap@cell.embeddings,
                                 as.character(seu@active.ident)))
    plot_df$UMAP_1 <- as.numeric(plot_df$UMAP_1)
    plot_df$UMAP_2 <- as.numeric(plot_df$UMAP_2)
    colnames(plot_df) <- c("UMAP_1","UMAP_2","CellType")


    # Plot diagnostic plot
    if(is.null(celltype_color)){
      celltype_color = colorRampPalette(brewer.pal(8, "Accent"))(length(levels(cell_label)))
      names(celltype_color) = levels(cell_label)
      # setting normal cells to black color
      celltype_color[gc] = "black"
    }
    if(is.null(celltype_cex)){ # point size
      celltype_cex=rep(2,length(levels(cell_label)))
    }
    if(is.null(celltype_pch)){ # point shape
      celltype_pch=1:length(levels(cell_label))
      celltype_pch[which(levels(cell_label)==gc)] = 0
    }

    if(plot_diagnostic==T){
      pdf(diagnostic_plot_path,width=20,height=16)
      cols=celltype_color#c("red","black","chartreuse4","orange")
      cex=celltype_cex # c(2,2,2,2)
      pchs = celltype_pch#c(1,3,2,1)
      par(mfrow=c(4,5))
    }
    # For each cluster
    pval = c()
    for (cluster in levels(cell_label)){
      if(cluster != "test"){ # skip normal cluster
        sel1 = which(cell_label==cluster)
        x1 = x[cell_label==cluster]
        n1= n[cell_label==cluster]
        # Binomial Test
        pnull=pjN[id]
        u1 = 1-pbinom(x1-1, n1, pnull)
        # FDR selection
        u1_order=order(u1)
        below_thresh = u1[u1_order]<FDR_alpha*c(1:length(u1))/length(u1)
        reject = sum(below_thresh)
        rejected=u1_order[which(below_thresh)]
        rejected_ind = sel1[rejected]

        # Add in UMAP information
        sig_cells = rep(0,dim(N)[2])
        sig_cells[rejected_ind] = 1 # 1 for significant, 0 for non-significant
        plot_df <- cbind(plot_df,sig_cells)
        pval = append(pval,u1)

        if(output==TRUE){
          # write the cell barcodes of the significant cells
          write.csv(data.frame(sig_cells=colnames(N)[rejected_ind]),
                  file=paste0(save_path,"/",voi[id],"_",
                              cluster,"_sig_cells.txt"),
                  row.names = FALSE,quote=FALSE)
          if(verbose==T){
            message(paste0("Significant cells of ", cluster," cluster saved!"))
          }
        }
        if(plot_diagnostic==T){
          # X vs N plot -- slope represents VAF
          plot(N[voi[id],], jitter(X[voi[id],], 0.25), col=cols[as.numeric(cell_label)],
            main=names(voi[id]), pch=pchs[as.numeric(cell_label)],
            cex=cex[as.numeric(cell_label)], xlab="Total Count", ylab="Alternative Allele Count")
          points(N[voi[id], rejected_ind], X[voi[id],rejected_ind],
              pch=18, cex=2, col=cols[as.numeric(unique(cell_label)[unique(cell_label) == cluster])])
          legend(x="topright",  col=c(cols[as.numeric(unique(cell_label))],
                                  cols[as.numeric(unique(cell_label)[unique(cell_label) == cluster])]),
              pch=c(pchs[as.numeric(unique(cell_label))],18),
              legend=paste(c(as.vector(unique(cell_label)),paste0(cluster," FDR 0.1"))),bg="white",bty="o")
          lines(N[voi[id],],N[voi[id],]*pnull, col="chartreuse4", lwd=2)
          text(x=max(N[voi[id],])*0.7,y=max(X[voi[id],])*(pnull+0.1),
            paste0("background freq=",round(pnull,6)),col="chartreuse4")
          grid()
          # X/N (VAF) vs N plot -- y represnts VAF
          plot(N[voi[id],], X[voi[id],]/(N[voi[id],]+0.1), xlab="Total Count",
            ylab="Alternative Allele Frequency", main=voi[id],
            col=cols[as.numeric(cell_label)], pch=pchs[as.numeric(cell_label)],
            cex=cex[as.numeric(cell_label)])
         points(N[voi[id], rejected_ind], X[voi[id],rejected_ind]/(N[voi[id],rejected_ind]+0.1),
              pch=18, cex=2, col=cols[as.numeric(unique(cell_label)[unique(cell_label) == cluster])])
          legend(x="topright",  col=c(cols[as.numeric(unique(cell_label))],
                                  cols[as.numeric(unique(cell_label)[unique(cell_label) == cluster])]),
             pch=c(pchs[as.numeric(unique(cell_label))],18),
             legend=paste(c(as.vector(unique(cell_label)),paste0(cluster," FDR 0.1"))),
             bg="white",bty="o")
          grid()
          abline(h=pnull, col="chartreuse4", lwd=2)
          text(x=max(N[voi[id],])*0.7,y=max(X[voi[id],]/(N[voi[id],]+0.1))*(pnull+0.1),
            paste0("background freq=",round(pnull,6)),col="chartreuse4")

          # Plot AF vs N profiles differently
          plot(N[voi[id],cell_label==cluster], X[voi[id],cell_label==cluster]/(N[voi[id],cell_label==cluster]+0.1),
            xlab="Total Count", ylab="Alternative Allele Frequency", main=paste0(cluster," AF vs Total Count"),
            col=cols[as.numeric(cell_label)][cell_label==cluster],
            pch=pchs[as.numeric(cell_label)][cell_label==cluster],
            cex=cex[as.numeric(cell_label)][cell_label==cluster],
            xlim=c(0,max(N[voi[id],])), ylim=c(0,1))
          points(N[voi[id], rejected_ind], X[voi[id],rejected_ind]/(N[voi[id],rejected_ind]+0.1),
              pch=18, cex=2, col=cols[as.numeric(unique(cell_label)[unique(cell_label) == cluster])])
          legend=paste(unique(cell_label))
          abline(h=pnull, col="chartreuse4", lwd=2)
          text(x=max(N[voi[id],])*0.7,y=pnull+0.1,
            paste0("background freq=",round(pnull,6)),col="chartreuse4")
          grid()
        }

       }else{ # normal cluster
         if(plot_diagnostic==T){
            plot(0,0, cex=0)
            plot(0,0, cex=0)
            plot(N[voi[id],cell_label==cluster], X[voi[id],cell_label==cluster]/(N[voi[id],cell_label==cluster]+0.1),
              xlab="Total Count", ylab="Alternative Allele Frequency", main=paste0(cluster," AF vs Total Count"),
              col=cols[as.numeric(cell_label)][cell_label==cluster],
              pch=pchs[as.numeric(cell_label)][cell_label==cluster],
              cex=cex[as.numeric(cell_label)][cell_label==cluster],
              xlim=c(0,max(N[voi[id],])), ylim=c(0,1))
            legend=paste(unique(cell_label))
            abline(h=pnull, col="chartreuse4", lwd=2)
            text(x=max(N[voi[id],])*0.7,y=pnull+0.1,
              paste0("background freq=",round(pnull,6)),col="chartreuse4")
          grid()

        # Diagnostic Histogram
        hist(af[cell_label==cluster & n>0], breaks=seq(0,1,0.05),
          xlab="allele frequency",
          col=cols[as.numeric(unique(cell_label)[unique(cell_label) == cluster])],
          main=cluster)
        # SQRT Histogram
        hist(sqrt(af[cell_label==cluster & n>0]), breaks=seq(0,1,0.05),
          xlab="SQRT allele frequency",
          col=cols[as.numeric(unique(cell_label)[unique(cell_label) == cluster])],
           main=cluster)
       }
      }
    }
    if(plot_diagnostic==T){
      dev.off()
    }
    if(verbose==T){
      message("Diagnostic Plot of the variant saved!")
    }
    pval_df= cbind(pval_df,pval)

    # Plotting significant cells' location in UMAP
    if(plot_umap==TRUE){
      colnames(plot_df) <- c(c("UMAP_1","UMAP_2","CellType"),levels(cell_label))
      g_list = list()
      k = 1
      for (cluster in levels(cell_label)){
        g <- ggplot(plot_df, aes(x = UMAP_1, y = UMAP_2)) +
        geom_point(size=1,col="gray") +
        geom_point(data=plot_df[plot_df[cluster] == 1,],aes(x = UMAP_1, y = UMAP_2),size=0.5,color="red")+
        theme(plot.title = element_text(hjust = 0.5)) +
        ggtitle(paste0("Signicant ",cluster, " Cells"))
        g_list[[k]] = g
        k=k+1
      }
      nrow=floor(sqrt(length(g_list)))
      ncol=ceiling(length(g_list)/nrow)
      pdf(umap_path,width=5*ncol,height=nrow*4)
      grid.arrange(grobs =g_list,nrow=nrow,top=voi[id])
      dev.off()
      if(verbose==T){
        message("UMAP of significant cells saved!")
      }
    }
  }
  colnames(pval_df) = voi
  return(pval_df)
}
