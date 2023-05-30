#!/usr/bin/Rscript

install.packages("BiocManager")
BiocManager::install(c('org.Hs.eg.db',"GOSemSim"))

# UniProt IDs
proteins=names(as.vector(read.csv("proteins.txt")[1,]))
fname="GO_similarity"

loadNamespace("Rcpp")

## Genome wide annotation for Human, primarily based on mapping using Entrez Gene identifiers
library("GOSemSim")
onthologies=c("BP","MF","CC")
hsGOs=c()
for (ont in onthologies) {
    print(ont)
    hsGOs <- c(hsGOs, godata('org.Hs.eg.db', keytype="UNIPROT", ont=ont, computeIC=TRUE))
}
names(hsGOs) <- onthologies

## Similarity as average of all three onthologies
similarity <- function(ont) {print(paste0(ont," similarity"));mgeneSim(genes=proteins, semData=hsGOs[[ont]], measure="Resnik", combine="BMA", verbose=T)}
for (ont in onthologies) {
    write.csv(x=similarity(ont), file=paste0(fname, "_", ont, ".csv"))
}

sim_MF <- read.csv(paste0(fname, "_MF.csv"), row.names=1)
sim_CC <- read.csv(paste0(fname, "_CC.csv"), row.names=1)
sim_BP <- read.csv(paste0(fname, "_BP.csv"), row.names=1)

reformat_mat <- function(sim, cols) {
    cols <- unique(cols)
    n <- length(cols)
    col_ls <- cols[cols %in% colnames(sim)]
    stopifnot(length(col_ls)>0)
    reformat_sim = data.frame(matrix(0, nrow=n, ncol=n), row.names=cols)
    colnames(reformat_sim) <- cols
    reformat_sim[col_ls, col_ls] = sim[col_ls, col_ls]
    reformat_sim
}

go_similarity = 1/3*(reformat_mat(sim_MF, proteins)+reformat_mat(sim_BP, proteins)+reformat_mat(sim_CC, proteins))
go_similarity[go_similarity==0] = NA
write.csv(x=go_similarity, file=paste0(fname, ".csv"))

system("rm -rf proteins.txt")
