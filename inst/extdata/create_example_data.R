# create example dataset from Recla et al. (2014) and Logan et al. (2013)
# See https://github.com/kbroman/qtl2data/tree/master/DO_Recla

library(qtl2geno)
json_file <- paste0("https://raw.githubusercontent.com/rqtl/",
                    "qtl2data/master/DO_Recla/recla.json")
recla <- read_cross2(json_file)

recla_pr <- calc_genoprob(recla, step=1, stepwidth="max", cores=0)
k <- calc_genetic_sim(recla_pr, use_grid_only=FALSE, cores=0)
y <- recla$pheno
X <- cbind(sex=as.numeric(recla$covar$Sex=="male"))
rownames(X) <- rownames(recla$covar)

# file to write a matrix to a CSV file
write_csv <-
    function(tab, file)
{
    tab <- cbind(ID=rownames(tab), tab)
    write.table(tab, file, row.names=FALSE, quote=FALSE,
                col.names=TRUE, sep=",")
}

# create temporary directory
dir <- tempdir()
dir <- file.path(dir, "Recla")
dir.create(dir)

# write files to that directory
write_csv(k, file=file.path(dir, "kinship.csv"))
write_csv(y, file=file.path(dir, "pheno.csv"))
write_csv(X, file=file.path(dir, "covar.csv"))

# write README
readme <- c("Example dataset from Recla et al. (2014) and Logan et al. (2013)",
            "See https://github.com/kbroman/qtl2data/tree/master/DO_Recla",
            "",
            "Contents (each with initial ID column):",
            "kinship.csv: kinship matrix",
            "pheno.csv:   phenotypes",
            "covar.csv:   sex covariate")
cat(paste0(readme, "\n"), file=file.path(dir, "ReadMe.txt"))

# create gzip file
wd <- getwd()
setwd(file.path(dir, ".."))
zip(file.path(wd, "recla.zip"), "Recla")
