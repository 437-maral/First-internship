corrected_refbase=function (beta = myNorm, arraytype = "450K") 
{
  message("[===========================]")
  message("[<<< ChAMP.REFBASE START >>>]")
  message("-----------------------------")
  projectWBC <- function(Y, coefWBC, contrastWBC = NULL, nonnegative = TRUE, 
                         lessThanOne = FALSE) {
    if (is.null(contrastWBC)) 
      Xmat = coefWBC
    else Xmat = coefWBC %*% t(contrastWBC)
    nCol = dim(Xmat)[2]
    nSubj = dim(Y)[2]
    mixCoef = matrix(0, nSubj, nCol)
    rownames(mixCoef) = colnames(Y)
    colnames(mixCoef) = colnames(Xmat)
    if (nonnegative) {
      if (lessThanOne) {
        Amat = cbind(rep(-1, nCol), diag(nCol))
        b0vec = c(-1, rep(0, nCol))
      }
      else {
        Amat = diag(nCol)
        b0vec = rep(0, nCol)
      }
      for (i in 1:nSubj) {
        obs = which(!is.na(Y[, i]))
        Dmat = t(Xmat[obs, ]) %*% Xmat[obs, ]
        mixCoef[i, ] = solve.QP(Dmat, t(Xmat[obs, ]) %*% 
                                  Y[obs, i], Amat, b0vec)$sol
      }
    }
    else {
      for (i in 1:nSubj) {
        obs = which(!is.na(Y[, i]))
        Dmat = t(Xmat[obs, ]) %*% Xmat[obs, ]
        mixCoef[i, ] = solve(Dmat, t(Xmat[obs, ]) %*% 
                               Y[obs, i])
      }
    }
    return(mixCoef)
  }
  message("<< Load projectWBC function success. >>")
  if (arraytype == "27K") {
    data(CellTypeMeans27K)
    DMRs27K <- intersect(rownames(CellTypeMeans27K), rownames(beta))
    cellFrac <- projectWBC(beta[DMRs27K, ], CellTypeMeans27K[DMRs27K, 
    ], lessThanOne = TRUE)
  }
  else {
    data(CellTypeMeans450K)
    DMRs450K <- intersect(rownames(CellTypeMeans450K), rownames(beta))
    cellFrac <- projectWBC(beta[DMRs450K, ], CellTypeMeans450K[DMRs450K, 
    ], lessThanOne = TRUE)
  }
  lm.o <- lm(t(beta) ~  cellFrac[,'Bcell'])
  tmp.m <- t(lm.o$res) + rowMeans(beta)
  tmp.m[tmp.m <= 0] <- min(tmp.m[which(tmp.m > 0)])
  tmp.m[tmp.m >= 1] <- max(tmp.m[which(tmp.m < 1)])
  message("All cell proportion influence except the one with least cell proportion get corrected.\n")
  message("[<<<< ChAMP.REFBASE END >>>>]")
  message("[===========================]")
  return(list(CorrectedBeta = tmp.m, CellFraction = cellFrac[,'Bcell']))
}
