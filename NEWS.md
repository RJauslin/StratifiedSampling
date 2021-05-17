# NEWS for Package "StratifiedSampling"
  
## StratitfiedSampling 0.1.9000
* add .9000 for dev version
* fix a findB when the number of columns is greater than the number of rows. 
* add calibRaking, harmonize, otmatch and bsmatch (a proposed method to do statistical matching)
* add vignettes to explain how to use matching.
* change README to add information about statistical matching.
  
## StratifiedSampling 0.1.0

### CRAN resubmission
* Please do not modify the user's global environment or the user's home filespace in your examples or vignettes by deleting objects:
rm(list = ls())
- remove rm(list = ls()) from R environment of cpp files. 


### CRAN resubmission
* moved https://arxiv.org/abs/arXiv:2101.05568 to https://arxiv.org/abs/2101.05568


### CRAN resubmission
* add \value to Rd file landingRM
* remove all rm(list = ls())
* remove all call to sampling:: and put sampling in suggest

### CRAN submission