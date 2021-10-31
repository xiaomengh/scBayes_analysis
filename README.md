scBayes manuscript additional analysis
======================================

This repository contains custom R code for additional analysis reported in the
scBayes manuscript. These analysis are not directly related to single cell
assignments. They are exploratory functional analyses after the cell identities
are assigned.

The repository contains two subdirectories, `BC_analysis` for the subclonal
ssGSEA analysis, and `CLL_analysis` for the single cell expression and
clonotype analyses.

To recreate the analyses results, navigate to the corresponding folder, and
source the R script from an R interpreter prompt or RStudio. Plots that are
parts of the figures included in the manuscript will be created, given you have
installed all the libraries used in the script. There is an additional step for
the CLL analyses before the R script can be sourced. The data file is too large
to be uploaded to github. We instead uploaded the file to amazon s3, and
provides a shell script to download the data. This script `fetch_data.sh` needs
to be run before the R script. In case the script fails, e.g. due to the
program `curl` not available on your system, you can also manually download the
data file using the following URL:

https://s3.us-west-2.amazonaws.com/public.marthlab.org/xiaomengh/Pt3_Bcells.rds
