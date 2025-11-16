## ----echo = FALSE-------------------------------------------------------------
  knitr::opts_chunk$set(
  eval=FALSE,
  results = "hide",
  collapse = TRUE, 
  comment = "#>",
  results = "asis"
)

## -----------------------------------------------------------------------------
#      # this will automatically download and install the dependencies.
#      install.packages("SNPannotator")

## -----------------------------------------------------------------------------
#  # first, load the library
#  library(SNPannotator)
#  
#  # save the configuration file to a specific folder
#  getConfigFile('/home/user/project1/postGWAS')

## -----------------------------------------------------------------------------
#  
#  library(SNPannotator)
#  
#  demo_annotation()
#  

## -----------------------------------------------------------------------------
#  
#  library(SNPannotator)
#  
#  getConfigFile('/home/user/config.ini')
#  

## -----------------------------------------------------------------------------
#  run_annotation('/home/user/analysis/config.ini')

## -----------------------------------------------------------------------------
#  findProxy("rs234")
#  
#  findProxy(c("rs234","rs678"))

## -----------------------------------------------------------------------------
#  output = findPairwiseLD(c('rs234','rs10244533'))
#  

## -----------------------------------------------------------------------------
#  findGenomicPos('rs121')
#  
#  findGenomicPos('7_24898067_A_G', type = "varid", file_path  = "output.xlsx")
#  

## -----------------------------------------------------------------------------
#  
#  findRSID(15, 79845218)
#  
#  findRSID(15, 79845218 , 79845238, file_path = "output.xlsx")
#  
#  findRSID(15, 80137560 ,build= "37")
#  

## -----------------------------------------------------------------------------
#  stringdb_annotation('BC_project',c('BRCA1','BRCA2'),limit=10)
#  

