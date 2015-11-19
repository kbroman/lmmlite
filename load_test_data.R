# data for testing lmmlite

if(!dir.exists("Rcache"))
    dir.create("Rcache")

file <- "Rcache/recla.RData" # to contain kinship matrix + phenotypes + covariates
if(file.exists(file)) {
    load(file)
} else {
    library(qtl2geno)
    json_file <- paste0("https://raw.githubusercontent.com/rqtl/",
                        "qtl2data/master/DO_Recla/recla.json")
    recla <- read_cross2(json_file)

    recla_pr <- calc_genoprob(recla, step=1, stepwidth="max", cores=0)
    k <- calc_genetic_sim(recla_pr, use_grid_only=FALSE, cores=0)
    y <- recla$pheno
    X <- cbind(intercept=1, sex=as.numeric(recla$covar$Sex=="male"))

    save(k, y, X, file=file)
}
