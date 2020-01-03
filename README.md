### Installation

**dmatch** relies on the following R packages: **Rcpp**, **RcppArmadillo**, **mvtnorm**, **RcppEigen**, **irlba**. All packagess are hosted on CRAN.

```r
install.packages("Rcpp")
install.packages("RcppArmadillo")
install.packages("mvtnorm")
install.packages("RcppEigen")
install.packages("irlba")
```

Install development version from GitHub:

```r
install.packages("devtools")
library(devtools)
devtools::install_github("qzhan321/dmatch")
```

### Quick start
For a quick start, visit our website to see a detailed tutorial with examples and the references for functions in **dmatch** illustrating its usage at:
https://qzhan321.github.io/dmatch/.

The complete analysis pipeline in the tutorial section is also included in this 
[R script](https://github.com/qzhan321/dmatch/blob/master/ExamplePipeline.R). 
Simply run:
```
R CMD BATCH ExamplePipeline.R
```
All the images and final corrected data will be generated in the current working directory. 

### Contact us
If you have any enquires, especially about performing dmatch integration on your data, please contact
**Qi Zhan** (UChicago) qizhan@uchicago.edu



