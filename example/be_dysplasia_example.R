library(Matrix)
library(scmtVT)

setwd("~/Documents/GitHub/scmtVT/")

# Step 1. load example data
seu <- readRDS("../data/a191_be/a191_be_final.rds")  # seurat object containing scRNA-Seq result and cluster annotations
counts <- read.table("../data/a191_be/a191_be_mt_counts.csv", sep=",", header=T, row.names = 1) # total varaint count matrix
af_mtx <- read.table("../data/a191_be/af_dm.csv",sep=",",header=T, row.names = 1) # alternate allele variant count matrix
vars_metadata <- read.table("../data/a191_be/191_be_allvariants.txt", sep="\t", header=T, row.names = 1) # metadata of the variants
# formatting data for ZIBB model 
N=as(as.matrix(counts[sapply(strsplit(rownames(af_mtx),"_"),"[[",1),]), "sparseMatrix") # total count coverage, variant x cell
rownames(N) <- rownames(af_mtx)
X =as(as.matrix(as.matrix(N)*af_mtx*0.01), "sparseMatrix")  # alternative allele count coverage, variant x cell
X = round(X)
# load cell labels and set background/normal cells
cell_label <- as.character(seu@active.ident)
cell_label[cell_label %in% c("FB","IM","VC")] = "Normal"
cell_label = as.factor(cell_label)

# load somatic variants of interest identified earlier
voi <- c("2692_G>A", "4037_G>A", "6360_G>A", "7074_G>A", "15153_G>A",
         "2813_T>C", "5055_T>C", "15305_T>C", "5215_T>C")

# Step 2. running expectation maximization to infer model parameters
gc = "Normal" # set background cluster name
X_sub = t(as.matrix(X[voi,cell_label == gc]))
N_sub = t(as.matrix(N[voi,cell_label == gc]))
# Expectation - Maximization to find the parameters from background/normal cells
set.seed(42)
param_fit=EM_ZIBB(X_sub,N_sub,alpha=1,beta=1,delta_gc=0.5,iterations=100,stop_diff=1e-05,gc="Normal",
                  save_path="example/output/zibb_fit/",verbose=T)
# (optional) create diagnostic plot
ZIBB_fit_diagnostic_plot(X_sub,N_sub,delta_gc=0.5,alpha=1,beta=1,delta_ij_hat=0.5,gc="Normal",
                                     delta_gc_list,LL_list,alpha_list,beta_list,save_path="example/output/zibb_fit/",
                                     width=10,height=10)
browseURL("example/output/zibb_fit/start_a_1_b_1_del_0.5_Normal.pdf")
# Step 3. Perform ZIBB statistical test for variants of interest
#  (optional) Re-load 100 iterations result
load("example/output/zibb_fit/start_a_1_b_1_del_0.5_Normal.RData")
# take parameter values from final iterations
alpha=alpha_list[length(alpha_list)]
beta=beta_list[length(beta_list)]
delta_gc=delta_gc_list[length(delta_gc_list)]
pval_df=ZIBB_test(X_sub,N_sub,alpha,beta,delta_ij_hat,gc="Normal",cell_label=cell_label,seu=seu,voi=voi,
                      output=T,plot_umap=T,save_path="example/output/zibb_fit/test_results/",verbose=F)
# not plotting
pval_df=ZIBB_test(X_sub,N_sub,alpha,beta,delta_ij_hat,gc="Normal",cell_label=cell_label,seu=seu,voi=voi,
                  output=F,plot_umap=T,save_path="example/output/zibb_fit/test_results/",verbose=F,plot_diagnostic =F)
# visualize one example
browseURL("example/output/zibb_fit/")

# Step 4. Calculate expected value of coverage-specific false discovery rate
FP_result = CoverageFP(save_path="example/output/zibb_fit/test_results/false_positives/",
                         X_sub,N_sub,alpha,beta,delta_ij_hat,bin_size = 20,voi,gc,FDR_alpha=0.1,cell_label,output=T)







