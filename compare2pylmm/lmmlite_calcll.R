# use lmmlite to calculate detailed log likelihood values

library(lmmlite)
data(recla)
n <- ncol(recla$pheno)
hsq_v <- seq(0, 1, len=101)
result <- NULL

cur <- 1
for(i in 1:ncol(recla$pheno)) {
    keep <- !is.na(recla$pheno[,i])
    y <- recla$pheno[keep,i,drop=FALSE]
    x <- recla$covar[keep,]
    k <- recla$kinship[keep,keep]
    n <- sum(keep)

    e <- eigen_rotation(k, y, x)
    res <- calcLL(hsq_v, e$Kva, e$y, e$X, reml=TRUE)
    result <- rbind(result,
                    data.frame(index=i, method="reml",
                               hsq=hsq_v,
                               loglik=res,
                               loglik_adj=res - n*(log(2*pi) + 1 - log(n))/2,
                               stringsAsFactors=FALSE))

    res <- calcLL(hsq_v, e$Kva, e$y, e$X, reml=FALSE)
    result <- rbind(result,
                    data.frame(index=i, method="ml",
                               hsq=hsq_v,
                               loglik=res,
                               loglik_adj=res - n*(log(2*pi) + 1 - log(n))/2,
                               stringsAsFactors=FALSE))
}

write.table(result, "lmmlite_llvals.csv", row.names=FALSE,
            col.names=TRUE, quote=FALSE, sep=",")
