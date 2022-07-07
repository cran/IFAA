## ---- include = FALSE---------------------------------------------------------
knitr::opts_chunk$set(
  collapse = TRUE,
  comment = "#>"
)

library(IFAA)

## ----eval=FALSE---------------------------------------------------------------
#  install.packages("IFAA", repos = "http://cran.us.r-project.org")

## ----eval=FALSE---------------------------------------------------------------
#  require(devtools)
#  devtools::install_github("gitlzg/IFAA")

## -----------------------------------------------------------------------------
library(IFAA)
suppressMessages(library(SummarizedExperiment))

## load the example microbiome data. This could be relative abundance or 
## absolute abundance data. If you have a csv or tsv file for the microbiome data, 
## you can use read.csv() function or read.table() function in R to read the 
## data file into R.
data(dataM)
dim(dataM)
dataM[1:5, 1:8]

## load the example covariates data. If you have a csv or tsv file for the 
## covariates data, you can use read.csv() function or read.table() function 
## in R to read the data file into R.
data(dataC)
dim(dataC)
dataC[1:3, ]

## -----------------------------------------------------------------------------
## Merge the microbiome data and covariate data by id to avoid unmatching observations. 
data_merged<-merge(dataM,dataC,by="id",all=FALSE)

## Seperate microbiome data and covariate data, drop id variable from microbiome data
dataM_sub<-data_merged[,colnames(dataM)[!colnames(dataM)%in%c("id")]]
dataC_sub<-data_merged[,colnames(dataC)]
 
## Create a SummarizedExperiment object 
test_dat<-SummarizedExperiment(assays=list(MicrobData=t(dataM_sub)), colData=dataC_sub)

## ---- eval=T------------------------------------------------------------------
results <- IFAA(experiment_dat = test_dat,
                testCov = c("v1"),
                ctrlCov = c("v2","v3"),
                sampleIDname = c("id"),
                fdrRate = 0.05)

## ----eval=T-------------------------------------------------------------------
summary_res<-results$full_result
sig_results<-subset(summary_res,sig_ind==T)
sig_results

## -----------------------------------------------------------------------------
## load the example microbiome data. This could be relative abundance or 
## absolute abundance data. If you have a csv or tsv file for the microbiome data, 
## you can use read.csv() function or read.table() function in R to read the 
## data file into R.
data(dataM)
dim(dataM)
dataM[1:5, 1:8]

## load the example covariates data. If you have a csv or tsv file for the 
## covariates data, you can use read.csv() function or read.table() function 
## in R to read the data file into R.
data(dataC)
dim(dataC)
dataC[1:3, ]

## -----------------------------------------------------------------------------
## load the example microbiome data. This could be relative abundance or 
## absolute abundance data. If you have a csv or tsv file for the microbiome data, 
## you can use read.csv() function or read.table() function in R to read the 
## data file into R.
data_merged<-merge(dataM,dataC,by="id",all=FALSE)

## load the covariates data. If you have a csv or tsv file for the covariates data, 
## you can use read.csv() function or read.table() function in R to read 
## the data file into R.
dataM_sub<-data_merged[,colnames(dataM)[!colnames(dataM)%in%c("id")]]
dataC_sub<-data_merged[,colnames(dataC)]
 
## Create a SummarizedExperiment object 
test_dat<-SummarizedExperiment(assays=list(MicrobData=t(dataM_sub)), colData=dataC_sub)

## ---- eval=T------------------------------------------------------------------
results <- MZILN(experiment_dat=test_dat,
                 refTaxa=c("rawCount11"),
                 allCov=c("v1","v2","v3"),
                 sampleIDname = c("id"),
                 fdrRate=0.05)

## -----------------------------------------------------------------------------
summary_res<-results$full_results

## -----------------------------------------------------------------------------
summary_res[summary_res$taxon=="rawCount18",,drop=FALSE]

## ----eval=T-------------------------------------------------------------------
subset(summary_res,sig_ind==T)

## -----------------------------------------------------------------------------
sessionInfo()

