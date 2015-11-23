# use lmmlite to fit LMM with the example data

library(lmmlite)
data(recla)
n <- ncol(recla$pheno)
result <- data.frame(index=rep(NA,2*n),
                     method=rep("", 2*n),
                     hsq=rep(NA,2*n),
                     intercept=rep(NA,2*n),
                     sex=rep(NA,2*n),
                     sigmasq=rep(NA,2*n),
                     loglik=rep(NA,2*n), stringsAsFactors=FALSE)

for(i in 1:ncol(recla$pheno)) {
    keep <- !is.na(recla$pheno[,i])
    y <- recla$pheno[keep,i,drop=FALSE]
    x <- recla$covar[keep,]
    k <- recla$kinship[keep,keep]

    e <- eigen_rotation(k, y, x)
    res <- fitLMM(e$Kva, e$y, e$X, reml=TRUE)
    result[2*i-1,1] <- i
    result[2*i-1,2] <- "reml"
    result[2*i-1,-(1:2)] <- c(res$hsq, res$beta, res$sigmasq, res$loglik)


    res <- fitLMM(e$Kva, e$y, e$X, reml=FALSE)
    result[2*i,1] <- i
    result[2*i,2] <- "ml"
    result[2*i,-(1:2)] <- c(res$hsq, res$beta, res$sigmasq, res$loglik)
}

write.table(result, "lmmlite_results.csv", row.names=FALSE,
            col.names=TRUE, quote=FALSE, sep=",")
