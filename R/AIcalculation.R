AIcalcu <- function(data,
                    ref,
                    Mprefix,
                    covsPrefix,
                    binPredInd,
                    contCovStd = FALSE) {
  results <- list()
  
  # get the original sample size
  nSub <- nrow(data)
  MVarNamLength <- nchar(Mprefix)
  
  # get taxa data
  micros <-
    lapply(substr(colnames(data), 1, MVarNamLength), function(x) {
      grep(Mprefix, x)
    })
  microPositions <- which(micros == 1)
  rm(micros)
  
  nTaxa <- length(microPositions)
  nNorm <- nTaxa - 1
  taxaNames <- colnames(data)[microPositions]
  rm(microPositions)
  
  # rearrange taxa names
  otherTaxaNames <- taxaNames[(taxaNames != ref)]
  taxaNames <- c(otherTaxaNames, ref)
  
  # get predictor data
  xVarNamLength <- nchar(covsPrefix)
  predics <-
    lapply(substr(colnames(data), 1, xVarNamLength), function(x) {
      grep(covsPrefix, x)
    })
  predPositions <- which(predics == 1)
  predNames <- colnames(data)[predPositions]
  nPredics <- length(predNames)
  rm(predics, predPositions)
  
  # taxa data
  w <- data[, taxaNames, drop = FALSE]
  
  # extract x data
  xData <- data[, predNames, drop = FALSE]
  
  if (contCovStd & (length(binPredInd) < nPredics)) {
    if (length(binPredInd) == 0) {
      xData <- scale(xData, center = FALSE)
    } else {
      xData[,-binPredInd] <-
        scale(xData[,-binPredInd, drop = FALSE], center = FALSE)
    }
  }
  
  rm(data, predNames)
  
  # transform data using log-ratio, creat Ai and Li
  l <- rep(NA, nSub)
  lLast <- rep(NA, nSub)
  logRatiow <- list()
  A <- list()
  for (i in seq_len(nSub)) {
    taxa.nonzero <- which(w[i,] != 0)
    if (length(taxa.nonzero) > 0) {
      lLast[i] <- max(taxa.nonzero)
    } else {
      lLast[i] <- 0
    }
    
    if (length(taxa.nonzero) > 0) {
      last.nonzero <- lLast[i]
      logwi <- as.numeric(log(w[i, taxa.nonzero]))
      l[i] <- length(logwi)
      if (l[i] > 1) {
        logRatiow[[i]] <- logwi[seq_len(l[i] - 1)] - logwi[l[i]]
        zero.m <- matrix(0, nrow = l[i] - 1, ncol = nNorm)
        if (last.nonzero == nTaxa) {
          zero.m[cbind(seq_len(l[i] - 1), taxa.nonzero[seq_len(l[i] - 1)])] <-
            1
        } else {
          zero.m[cbind(seq_len(l[i] - 1), taxa.nonzero[seq_len(l[i] - 1)])] <-
            1
          zero.m[, taxa.nonzero[l[i]]] <- -1
        }
        
        A[[i]] <- MatrixExtra::as.csc.matrix(zero.m)
        rm(zero.m)
      } else {
        logRatiow[[i]] <- NA
        A[[i]] <- NA
      }
    } else {
      l[i] <- 0
      logRatiow[[i]] <- NA
      A[[i]] <- NA
    }
  }
  
  # obtain the list of samples whose have at least 2 non-zero taxa
  twoList <- which(l > 1)
  lengthTwoList <- length(twoList)
  
  rm(w)
  
  results$xData <- xData
  rm(xData)
  
  results$logRatiow <- logRatiow
  rm(logRatiow)
  results$A <- A
  rm(A)
  results$twoList <- twoList
  rm(twoList)
  results$taxaNames <- taxaNames
  rm(taxaNames)
  results$lengthTwoList <- lengthTwoList
  results$lLast <- lLast
  results$l <- l
  results$nTaxa <- nTaxa
  results$nNorm <- nNorm
  results$nSub <- nSub
  results$nPredics <- nPredics
  return(results)
}
