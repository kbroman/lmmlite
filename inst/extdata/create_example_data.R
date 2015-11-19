# Script to create the example dataset
# The dataset is from Recla et al. (2014) and Logan et al. (2013)
# See https://github.com/kbroman/qtl2data/tree/master/DO_Recla

library(qtl2geno)
json_file <- paste0("https://raw.githubusercontent.com/rqtl/",
                    "qtl2data/master/DO_Recla/recla.json")
recla <- read_cross2(json_file)

recla_pr <- calc_genoprob(recla, step=1, stepwidth="max", cores=0)
k <- calc_genetic_sim(recla_pr, use_grid_only=FALSE, cores=0)
y <- recla$pheno
X <- cbind(intercept=1, sex=as.numeric(recla$covar$Sex=="male"))
rownames(X) <- rownames(recla$covar)

recla <- list(kinship=k, pheno=y, covar=X)
save(recla, file="recla.RData")
