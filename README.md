# scmtVT

This GitHub contains the beta version of the Single Cell Mitochondrial Variant Test (scmtVT) package.<br/>
It implements a zero-inflated beta-binomial (ZIBB) test to identify cells that are significant for a given variant of interest in scRNA-seq with mitochondrial enrichment data.

### Installation
Dependency of this package is VGAM. 
Please install the dependencies first. The following commands could be used:
``` r
install.packages("VGAM")
```

To install the package, please use the following commands:
``` r
install.packages("devtools")
devtools::install_github("jiazhen-rong/scmtVT") # install
library(scmtVT) # load
```
or directly copy from git:
``` linux
git clone https://github.com/jiazhen-rong/scmtVT.git
```
### Tutorial

An example of how to run the ZIBB test is shown below:

- [Tutorial for the ZIBB test on a dysplastic Barrett's esophagus sample](https://github.com/jiazhen-rong/scmtVT/blob/master/example/)

### Citations
If you used the package in your research, please cite:

*Clonal cell states link Barrett’s esophagus and esophageal adenocarcinoma
Rodrigo A. Gier, Raúl A. Reyes Hueros, Jiazhen Rong, Maureen DeMarshall, Tatiana A. Karakasheva, Amanda B. Muir, Gary W. Falk, Nancy R. Zhang, Sydney M. Shaffer
bioRxiv 2023.01.26.525564; doi: https://doi.org/10.1101/2023.01.26.525564*

### License
Apache License 2.0
