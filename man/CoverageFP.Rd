% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/CoverageFDR.R
\name{CoverageFP}
\alias{CoverageFP}
\title{This function is for Coverage-specific False Positives Estimation,
calculating the number of expected false positives in the significant cells.}
\usage{
CoverageFP(
  X_sub,
  N_sub,
  alpha,
  beta,
  delta_ij_hat,
  bin_size = 20,
  voi,
  gc,
  FDR_alpha = 0.1,
  cell_label,
  output = T,
  save_path = NULL
)
}
\arguments{
\item{X_sub}{alternate allele count matrix.}

\item{N_sub}{total variant count matrix.}

\item{alpha}{parameter of beta prior, estimated from EM.}

\item{beta}{parameter of beta prior, estimated from EM.}

\item{delta_ij_hat}{parameter representing non-zero distribution component, estimated from EM.}

\item{bin_size}{The bin size of coverage. Default is 20 counts per bin.}

\item{voi}{A list of Variant of Interest after basic filtering.}

\item{gc}{name of background cluster, used for saving the file.}

\item{FDR_alpha}{A numeric value representing FDR threshold.}

\item{cell_label}{A vector containing cell identities based on scRNA-Seq and standard Seurat annotation}

\item{output}{TRUE or FALSE, representing whether to write significant cells of each cluster.}

\item{save_path}{path to save the diagnostic plots. Default setting will store in "zibb_fit" folder under current working directory.}

\item{seu}{A seurat object generated from standard scRNA-Seq processing.}

\item{plot_umap}{TRUE or FALSE, representing whether to plot significant cells of each cluster on UMAP.}

\item{plot_diagnostic}{TRUE or FALSE, representing whether to plot diagnostic plot.}

\item{verbose}{TRUE or FALSE, representing whether to print out detailed log information.}
}
\value{
A list containing:
\itemize{
\item sig_cell_df Dataframe containing number of significant cells in each variant.
\item FP_df Dataframe containing number of expected false positives in each variant.
}
}
\description{
This function is for Coverage-specific False Positives Estimation,
calculating the number of expected false positives in the significant cells.
}
\examples{
pval_df=ZIBB_test(X_sub,N_sub,alpha,beta,delta_ij_hat,gc="Normal",cell_label=cell_label,seu=seu,
          output=T,plot_umap=T,save_path="example/output/zibb_fit/plots/")

}
