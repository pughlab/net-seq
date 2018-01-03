.libPaths(c(.libPaths(), "/mnt/work1/users/pughlab/bin/r_lib/library/3.1"))
#################################################
#######
#######       Sample-GTEx 
#######       Usage: Rscript sampleToGtexCdf.R /path/to/sampleTpm.Rdata /path/to/goi_list.txt outdir
#######
#######   Author: Rene Quevedo
#######   Date: July 27, 2016
#######   Notes: FPKM processed by Cindy
#######        : FPKM and TPM matrices assembled by Rene
#################################################
args <- commandArgs(TRUE)

load("/mnt/work1/users/pughlab/projects/NET-SEQ/loh_signature/wes_validation-cohort/data/loh-het_regions/NETseq_loh-het-genes.Rdata")
load("/mnt/work1/users/pughlab/projects/NET-SEQ/rna_seq/Cufflinks/star2.4.2a/data_matrix/tpm/rdata/NETSEQ_on-treatment.annotation.Rdata")
load("/mnt/work1/users/pughlab/projects/NET-SEQ/rna_seq/Cufflinks/star2.4.2a/data_matrix/tpm/rdata/NETSEQ_pre-treatment.tpm.Rdata")
pre.tpm.mat <- tpm.mat
load("/mnt/work1/users/pughlab/projects/NET-SEQ/rna_seq/Cufflinks/star2.4.2a/data_matrix/tpm/rdata/NETSEQ_on-treatment.tpm.Rdata")
sample.mat <- cbind(tpm.mat, pre.tpm.mat)
# [1] "global.anno.df" "het.genes"      "loh.genes"      "pre.tpm.mat"   "tpm.mat"
goi.v <- loh.genes$genes
goi.v <- factor(c(as.character(goi.v), 'AURKA'))

# Temporary GTEx Source Files
source('/mnt/work1/users/pughlab/bin/gtex_tools/src/tpmPreprocess.R')
source('/mnt/work1/users/pughlab/bin/gtex_tools/src/visualizeTpm.R')
load('/mnt/work1/users/pughlab/bin/gtex_tools/src/annoDf.RData')

#GTEx Files on Mordor
all.tissue.id <- '/mnt/work1/users/pughlab/external_data/GTEx/FPKMmatrix/tpm/rdata/allTissue.txt'
tpm.dir <- '/mnt/work1/users/pughlab/external_data/GTEx/FPKMmatrix/tpm/rdata'
master.annotation.file <- '/mnt/work1/users/pughlab/external_data/GTEx/Annotations/gtex.sra.annot.histFormatted.RData'
all.tissues <- read.table(all.tissue.id, header=FALSE, stringsAsFactors=FALSE, check.names=FALSE)
load(master.annotation.file)

# Compare each sample in sample.tpm to each tissue sample in tpm.tissue.list
#load(sample.tpm)
#sample.mat <- tpm.mat
#goi.v <- read.table(goi.txt, header=FALSE, check.names=FALSE, stringsAsFactors=FALSE)$V1
# Creates a list of TPM matrices given a list of tissues of interest and extracts only the given gene
tpm.tissue.list <- splitTissueTpm(all.tissues=all.tissues$V1)

# Convert the genes of interest to tracking_ids and subsets matrices based on this
#goi.v <- c("DNMT1", "CCNE1", "CCNE2", "SOX2", "CDKN1B", "CDKN1A", "MEN1", "CCND1", "CDK2", "CDK4", "RB1", "CTNNB1", "DAXX", "ATRX")
goi.tracking_id.v <- sapply(goi.v, function(gene.name) getId(gene.name, id.type='hugo', global.anno.df=global.anno.df))
names(goi.tracking_id.v) <- goi.v




ecdf.list <- list()
ci95.list <- list()
quant.list <- list()
net.list <- list()
for(each.goi.num in c(1:length(goi.tracking_id.v))){
  goi.hugo <- names(goi.tracking_id.v)[each.goi.num]
  each.goi <- goi.tracking_id.v[each.goi.num]
  
  print(paste("Processing GOI: ", goi.hugo, sep=""))
  # Isolate GTEx sample ecdf for each GOI
  tpm.tissue.subset.ecdf <- lapply(tpm.tissue.list, function(tissue.tpm.mat){
    tissue.tpm.mat <- tissue.tpm.mat[which(rownames(tissue.tpm.mat) %in% each.goi),,drop=FALSE]
    tissue.tpm.sub.mat <- tissue.tpm.mat[match(each.goi, rownames(tissue.tpm.mat)),,drop=FALSE]
    
    tissue.tpm.sub.95ci <- get95Ci(tissue.tpm.sub.mat[1,])  #Tissue and gene specific 95% CI
    tissue.tpm.sub.ecdf <- tryCatch({ ecdf(tissue.tpm.sub.mat[1,]) # Tissue and gene specific ECDF
    }, error=function(e){ NA })
    
    return(list(ci95=tissue.tpm.sub.95ci,
                ecdf=tissue.tpm.sub.ecdf))
  })
  
  # Isolate the tumor sample expression percentile in comparison to GTEx Tissue/gene ECDF
  sample.subset.mat <- sample.mat[which(rownames(sample.mat) %in% each.goi),,drop=FALSE]
  sample.subset.mat <- sample.subset.mat[match(each.goi, rownames(sample.subset.mat)),,drop=FALSE]
  tissue.ecdf.list <- lapply(tpm.tissue.subset.ecdf, function(each.ecdf) {
    ecdf.sample.mat <- tryCatch({each.ecdf[['ecdf']](sample.subset.mat)
    }, error=function(e){matrix(rep(NA, ncol(sample.subset.mat)), nrow=1)})
    names(ecdf.sample.mat) <- colnames(sample.subset.mat)
    return(ecdf.sample.mat)
  })
  tissue.ecdf.df <- do.call("rbind", tissue.ecdf.list)
  if(!all(is.na(tissue.ecdf.df))){
    hi.quant.ecdf <- apply(tissue.ecdf.df, 2, function(x) quantile(x, 0.95))
    pnet.idx <- unlist(sapply(c("001", "003", "008", "009"), function(x) grep(x, colnames(tissue.ecdf.df)) ))
    ginet.idx <- seq(1:ncol(tissue.ecdf.df))[-pnet.idx]
    quant.list[[goi.hugo]] <- list("pnet"=hi.quant.ecdf[pnet.idx],
                                   "ginet"=hi.quant.ecdf[ginet.idx])
    
    mean.gtex.expr <- sapply(tpm.tissue.subset.ecdf, function(x) x$ci95[1])
    net.list[[goi.hugo]] <- c("med_PNET"=median(sample.subset.mat[pnet.idx]),
                              "mean_PNET"=mean(sample.subset.mat[pnet.idx]),
                              "med_GINET"=median(sample.subset.mat[ginet.idx]),
                              "mean_GINET"=mean(sample.subset.mat[ginet.idx]),
                              "gtex"=mean(mean.gtex.expr))
  }
  
  
  
  # Append the list containing all the Percentile values and 95% CI values
  # ecdf.list[[goi.hugo]] <- tissue.ecdf.list[order(sapply(tissue.ecdf.list, median))]
  # ci95.list[[goi.hugo]] <- lapply(tpm.tissue.subset.ecdf[order(sapply(tissue.ecdf.list, median))],
  #                                 function(x) x[['ci95']])
}
pnet.quant.df <- as.data.frame(do.call("rbind", lapply(quant.list, function(x) x[['pnet']])))
ginet.quant.df <- as.data.frame(do.call("rbind", lapply(quant.list, function(x) x[['ginet']])))
pnet.quant.df$meanQ <- apply(pnet.quant.df, 1, mean)
ginet.quant.df$meanQ <- apply(ginet.quant.df, 1, mean)

low.expr.idx <- pnet.quant.df$meanQ < 0.15
pnet.lowExpr <- pnet.quant.df[low.expr.idx,]
ginet.lowExpr <- ginet.quant.df[low.expr.idx,]

do.call("rbind", net.list)
#
#aurka.pancreas.ecdf <- ecdf(gtex.pancreas.mat[grep("ENSG00000087586", rownames(gtex.pancreas.mat)),])
#pnet.aurka.ecdf <- aurka.pancreas.ecdf(tpm.mat[1768,])
# pnet.ecdf.vals <- pnet.aurka.ecdf[pnet.tpm.idx]
# names(pnet.ecdf.vals) <- colnames(tpm.mat)[pnet.tpm.idx]
# ginet.ecdf.vals <- pnet.aurka.ecdf[ginet.tpm.idx]
# names(ginet.ecdf.vals) <- colnames(tpm.mat)[ginet.tpm.idx]

save(quant.list, )
write.table(pnet.quant.df[low.expr.idx,], file="pnet2gtex_lowExpr.tsv",
            sep="\t", col.names = TRUE, row.names = TRUE, quote=FALSE)
write.table(ginet.quant.df[low.expr.idx,], file="ginet2gtex_lowExpr.tsv",
            sep="\t", col.names = TRUE, row.names = TRUE, quote=FALSE)