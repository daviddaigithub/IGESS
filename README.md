IGESS
=======

IGESS is a statistical approach to integrating individual level genotype data and summary statistics in Genome Wide Association Studies. 'IGESS' R package provides computationally efficient and user friendly interface to fit and evaluate the IGESS model. It accepts both the R-type data  and binary plink files.

Usage
=======

The file IGESS_package.pdf in inst/doc shows several examplex of how to use IGESS package. 

Development 
=======
This R package is developed by Mingwei Dai and Can Yang, and maintained by Can Yang <eeyangc@gmail.com>.

Installation
=======
To install the development version of IGESS, it's easiest to use the 'devtools' package. Note that IGESS depends on the 'Rcpp' and 'RcppArmadillo' package, which also requires appropriate setting of Rtools and Xcode for Windows and Mac OS/X, respectively.

install.packages("devtools")  
library(devtools)  
install_github("daviddaigithub/IGESS")  

References
=======
M. Dai, J. Ming, M., Cai, J. Liu, C. Yang, X. Wan, and Z. Xu. IGESS: A statistical approach to integrating individual level genotype data and summary statistics in genome wide association studies. Bioinformatics. 2017 Sep 15;33(18):2882-2889. doi: 10.1093/bioinformatics/btx314

=======
# IGESS
IGESS(Gentic Analysis integrating individual level data and summary statistics)

