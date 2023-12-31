% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/ZIBBTest.R
\name{ZIBB_test}
\alias{ZIBB_test}
\title{This function contains corresponding ZIBB test for variant of interest
and creates corresponding diagnostic plots.
The significant cells for each variant of interest will be labeled in the diagnostic plot and the UMAP.}
\usage{
ZIBB_test(
  X_sub,
  N_sub,
  alpha,
  beta,
  delta_ij_hat,
  FDR_alpha = 0.1,
  gc = "Normal",
  cell_label = NULL,
  seu = NULL,
  voi = NULL,
  output = T,
  plot_umap = T,
  save_path = NULL,
  plot_diagnostic = T,
  verbose = F
)
}
\arguments{
\item{X_sub}{alternate allele count matrix.}

\item{N_sub}{total variant count matrix.}

\item{alpha}{parameter of beta prior, estimated from EM.}

\item{beta}{parameter of beta prior, estimated from EM.}

\item{delta_ij_hat}{parameter representing non-zero distribution component, estimated from EM.}

\item{FDR_alpha}{A numeric value representing FDR threshold.}

\item{gc}{name of background cluster, used for saving the file.}

\item{cell_label}{A vector containing cell identities based on scRNA-Seq and standard Seurat annotation}

\item{seu}{A seurat object generated from standard scRNA-Seq processing.}

\item{voi}{A list of Variant of Interest after basic filtering.}

\item{output}{TRUE or FALSE, representing whether to write significant cells of each cluster.}

\item{plot_umap}{TRUE or FALSE, representing whether to plot significant cells of each cluster on UMAP.}

\item{save_path}{path to save the diagnostic plots. Default setting will store in "zibb_fit" folder under current working directory.}

\item{plot_diagnostic}{TRUE or FALSE, representing whether to plot diagnostic plot.}

\item{verbose}{TRUE or FALSE, representing whether to print out detailed log information.}
}
\value{
A matrix containing ZIBB test's p-values of each cell i for each given variant of interest j (cell x variant)
}
\description{
This function contains corresponding ZIBB test for variant of interest
and creates corresponding diagnostic plots.
The significant cells for each variant of interest will be labeled in the diagnostic plot and the UMAP.
}
\examples{
pval_df=ZIBB_test(X_sub,N_sub,alpha,beta,delta_ij_hat,gc="Normal",cell_label=cell_label,seu=seu,
          output=T,plot_umap=T,save_path="example/output/zibb_fit/plots/")

}
