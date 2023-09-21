# scmtVT

This Github contains the beta version package for Single Cell Mitochondrial Variant Test (scmtVT).<br/>
It implements a zero-inflated beta-binomail test for identifying the significant carriers of a given variant of interest.

### Installation

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

An example of identifying significant carrier through ZIBB test is shown below:

- [Tutorial of ZIBB test on Dysplasia and Barrett's Esaphugus example](https://github.com/seasoncloud/Clonalscope/tree/main/samples/P5931/scRNA)

### Citations
If you used the package for formal studies, please cite the following paper:

Clonal cell states link Barrett’s esophagus and esophageal adenocarcinoma
Rodrigo A. Gier, Raúl A. Reyes Hueros, Jiazhen Rong, Maureen DeMarshall, Tatiana A. Karakasheva, Amanda B. Muir, Gary W. Falk, Nancy R. Zhang, Sydney M. Shaffer
bioRxiv 2023.01.26.525564; doi: https://doi.org/10.1101/2023.01.26.525564

### Lisences
Apache 2.0 Liscence
