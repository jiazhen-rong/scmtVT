# scmtVT

### Overview
This GitHub contains the beta version of the Single Cell Mitochondrial Variant Test (scmtVT) package.<br/> It implements a zero-inflated beta-binomial (ZIBB) test to identify cells that are significant for a given variant of interest in scRNA-seq with mitochondrial enrichment data.

All the results in the paper is based on the original version v1.0.0.

### Updates 2025-08-20

The diagnostic plot of each patient is now available in the [results](https://github.com/jiazhen-rong/scmtVT/blob/master/results/) folder.

### Systems Requirement
#### Hardware requirements
`scmtVT` requires a standard computer with enough RAM and memory.
#### Software requirements
##### OS requirement
The software was tested on:
   - MacOS Sequoia (15.6.1)
   - Rocky Linux 8.10 (Green Obsidian)
##### R dependencies
The package depends on the following R packages:
```
RColorBrewer, VGAM, ggplot2, gridExtra, cli,         
glue, grDevices, grid, gtable, isoban,   
lifecycle, MASS, mgcv, rlang, scales,      
stats, tibble, vctrs, withr, graphics,    
utils, methods, stats4, splines, nlme,        
Matrix, farver, labeling, munsell, R6,          
viridisLite, fansi, magrittr, pillar, pkgconfig,   
lattice, colorspace, utf8
```

### Installation
The installation of the package shall take < 5min. If question was encountered during installation, feel free to post in Issues.

Please install VGAM first. 

The following commands could be used:

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

-   [Tutorial for the ZIBB test on a dysplastic Barrett's esophagus sample](https://github.com/jiazhen-rong/scmtVT/blob/master/example/). This demo will take ~ 20min. 

### Citations

If you used the package in your research, please cite:

*Clonal cell states link Barrett's esophagus and esophageal adenocarcinoma Rodrigo A. Gier, Raúl A. Reyes Hueros, Jiazhen Rong, Maureen DeMarshall, Tatiana A. Karakasheva, Amanda B. Muir, Gary W. Falk, Nancy R. Zhang, Sydney M. Shaffer bioRxiv 2023.01.26.525564; doi: <https://doi.org/10.1101/2023.01.26.525564>*

### License

Apache License 2.0
