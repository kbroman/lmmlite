library(lmmlite)
data(recla)
for(i in c("kinship", "pheno", "covar"))
    write.table(recla[[i]], paste0(i, ".csv"), sep=",", quote=FALSE,
                row.names=FALSE, col.names=FALSE)
