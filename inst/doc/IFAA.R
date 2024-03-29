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

set.seed(1)

## create an ID variable for the example data
ID <- seq_len(20)

## generate three covariates x1, x2, and x3, with x2 binary
x1 <- rnorm(20)
x2 <- rbinom(20, 1, 0.5)
x3 <- rnorm(20)
dataC <- data.frame(cbind(ID, x1, x2, x3))

## Coefficients for x1, x2, and x3 among 10 taxa.
beta_1 <- c(0.1, rep(0, 9))
beta_2 <- c(0, 0.2, rep(0, 8))
beta_3 <- rnorm(10)
beta_mat <- cbind(beta_1, beta_2, beta_3)

## Generate absolute abundance for 10 taxa in ecosystem.
dataM_eco <- floor(exp(10 + as.matrix(dataC[, -1]) %*% t(beta_mat) + rnorm(200, sd = 0.05)))

## Generate sequence depth and generate observed abundance
Ci <- runif(20, 0.01, 0.05)
dataM <- floor(apply(dataM_eco, 2, function(x) x * Ci))
colnames(dataM) <- paste0("rawCount", 1:10)

## Randomly introduce 0 to make 25% sparsity level.
dataM[sample(seq_len(length(dataM)), length(dataM) / 4)] <- 0

dataM <- data.frame(cbind(ID, dataM))

## The following steps are to create a SummarizedExperiment object.
## If you already have a SummarizedExperiment format data, you can
## ignore the following steps and directly feed it to the IFAA function.

## Merge two dataset by ID variable
data_merged <- merge(dataM, dataC, by = "ID", all = FALSE)

## Seperate microbiome data and covariate data, drop ID variable from microbiome data
dataM_sub <- data_merged[, colnames(dataM)[!colnames(dataM) %in% c("ID")]]
dataC_sub <- data_merged[, colnames(dataC)]

## Create SummarizedExperiment object
test_dat <- SummarizedExperiment::SummarizedExperiment(
  assays = list(MicrobData = t(dataM_sub)), colData = dataC_sub
)

## ---- eval=TRUE---------------------------------------------------------------
set.seed(123) # For full reproducibility
results <- IFAA(
  experiment_dat = test_dat,
  testCov = c("x1", "x2"),
  ctrlCov = c("x3"),
  sampleIDname = "ID",
  fdrRate = 0.05,
  nRef = 2,
  paraJobs = 2
)

## ----eval=TRUE----------------------------------------------------------------
summary_res <- results$full_results
sig_results <- subset(summary_res, sig_ind == TRUE)
sig_results

## ---- eval=FALSE--------------------------------------------------------------
#  set.seed(123) # For full reproducibility
#  
#  results <- IFAA(
#    experiment_dat = test_dat,
#    microbVar = c("rawCount1", "rawCount2", "rawCount3"),
#    testCov = c("x1", "x2"),
#    ctrlCov = c("x3"),
#    sampleIDname = "ID",
#    fdrRate = 0.05,
#    nRef = 2,
#    paraJobs = 2
#  )

## ---- eval=TRUE---------------------------------------------------------------
set.seed(123) # For full reproducibility

resultsRatio <- MZILN(
  experiment_dat = test_dat,
  microbVar = c("rawCount5","rawCount8"),
  refTaxa = c("rawCount10"),
  allCov = c("x1", "x2", "x3"),
  sampleIDname = "ID",
  fdrRate = 0.05,
  paraJobs = 2
)

## -----------------------------------------------------------------------------
summary_res_ratio <- resultsRatio$full_results
summary_res_ratio

## ---- eval=T------------------------------------------------------------------
set.seed(123) # For full reproducibility
resultsAllRatio <- MZILN(
  experiment_dat = test_dat,
  refTaxa = c("rawCount10"),
  allCov = c("x1", "x2", "x3"),
  sampleIDname = "ID",
  fdrRate = 0.05,
  paraJobs = 2
)

## -----------------------------------------------------------------------------
summary_res_ratio <- resultsAllRatio$full_results
summary_res_ratio[summary_res_ratio$taxon == "rawCount5", , drop = FALSE]

## -----------------------------------------------------------------------------
sessionInfo()

