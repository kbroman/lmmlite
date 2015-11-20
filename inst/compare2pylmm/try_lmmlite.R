# use lmmlite to fit LMM with the example data

library(lmmlite)
data(recla)
n <- ncol(recla$pheno)
result <- data.frame(index=rep(NA,2*n),
                     method=rep("", 2*n),
                     hsq=rep(NA,2*n),
                     intercept=rep(NA,2*n),
                     sex=rep(NA,2*n),
                     sigsq=rep(NA,2*n),
                     loglik=rep(NA,2*n), stringsAsFactors=FALSE)

cur <- 1
for(i in 1:ncol(recla$pheno)) {
    keep <- !is.na(recla$pheno[,i])
    y <- recla$pheno[keep,i,drop=FALSE]
    x <- recla$covar[keep,]
    k <- recla$kinship[keep,keep]

    e <- eigen_rotation(k, y, x)
    res <- fitLMM(e$Kva, e$y, e$X, reml=TRUE)
    result[cur,1] <- i
    result[cur,2] <- "reml"
    result[cur,-(1:2)] <- c(res$hsq, res$beta, res$sigsq, res$loglik)

    cur <- cur + 1
    e <- eigen_rotation(k, y, x)

    res <- fitLMM(e$Kva, e$y, e$X, reml=FALSE)
    result[cur,1] <- i
    result[cur,2] <- "ml"
    result[cur,-(1:2)] <- c(res$hsq, res$beta, res$sigsq, res$loglik)

    cur <- cur + 1
}

write.table(result, "lmmlite_results.csv", row.names=FALSE,
            col.names=TRUE, quote=FALSE, sep=",")
