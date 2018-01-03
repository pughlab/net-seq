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
sample.tpm <- args[1]   # Rdata dataframe of tpm, sample x tracking_id
goi.txt <- args[2]      # Genes of interest
outdir <- args[3]      # Genes of interest

# Temporary GTEx Source Files
source('~/git/net-seq/gtex_analysis/src/tpmPreprocess.R')
source('~/git/net-seq/gtex_analysis/src/visualizeTpm.R')
load('~/git/net-seq/gtex_analysis/src/annoDf.RData')

#GTEx Files on Mordor
all.tissue.id <- '/mnt/work1/users/pughlab/external_data/GTEx/FPKMmatrix/tpm/rdata/allTissue.txt'
tpm.dir <- '/mnt/work1/users/pughlab/external_data/GTEx/FPKMmatrix/tpm/rdata'
master.annotation.file <- '/mnt/work1/users/pughlab/external_data/GTEx/Annotations/gtex.sra.annot.histFormatted.RData'
all.tissues <- read.table(all.tissue.id, header=FALSE, stringsAsFactors=FALSE, check.names=FALSE)
load(master.annotation.file)

# Compare each sample in sample.tpm to each tissue sample in tpm.tissue.list
load(sample.tpm)
sample.mat <- tpm.mat
goi.v <- read.table(goi.txt, header=FALSE, check.names=FALSE, stringsAsFactors=FALSE)$V1
# Creates a list of TPM matrices given a list of tissues of interest and extracts only the given gene
tpm.tissue.list <- splitTissueTpm(all.tissues=all.tissues$V1)

# Convert the genes of interest to tracking_ids and subsets matrices based on this
#goi.v <- c("DNMT1", "CCNE1", "CCNE2", "SOX2", "CDKN1B", "CDKN1A", "MEN1", "CCND1", "CDK2", "CDK4", "RB1", "CTNNB1", "DAXX", "ATRX")
goi.tracking_id.v <- sapply(goi.v, function(gene.name) getId(gene.name, id.type='hugo', global.anno.df=global.anno.df))
names(goi.tracking_id.v) <- goi.v



ecdf.list <- list()
ci95.list <- list()
for(each.goi.num in c(1:length(goi.tracking_id.v))){
  goi.hugo <- names(goi.tracking_id.v)[each.goi.num]
  each.goi <- goi.tracking_id.v[each.goi.num]
  
  print(paste("Processing GOI: ", goi.hugo, sep=""))
  # Isolate GTEx sample ecdf for each GOI
  tpm.tissue.subset.ecdf <- lapply(tpm.tissue.list, function(tissue.tpm.mat){
    tissue.tpm.mat <- tissue.tpm.mat[which(rownames(tissue.tpm.mat) %in% each.goi),,drop=FALSE]
    tissue.tpm.sub.mat <- tissue.tpm.mat[match(each.goi, rownames(tissue.tpm.mat)),,drop=FALSE]
    
    tissue.tpm.sub.95ci <- get95Ci(tissue.tpm.sub.mat[1,])  #Tissue and gene specific 95% CI
    tissue.tpm.sub.ecdf <- ecdf(tissue.tpm.sub.mat[1,]) # Tissue and gene specific ECDF
    return(list(ci95=tissue.tpm.sub.95ci,
                ecdf=tissue.tpm.sub.ecdf))
  })
  
  # Isolate the tumor sample expression percentile in comparison to GTEx Tissue/gene ECDF
  sample.subset.mat <- sample.mat[which(rownames(sample.mat) %in% each.goi),,drop=FALSE]
  sample.subset.mat <- sample.subset.mat[match(each.goi, rownames(sample.subset.mat)),,drop=FALSE]
  tissue.ecdf.list <- lapply(tpm.tissue.subset.ecdf, function(each.ecdf) {
    ecdf.sample.mat <- each.ecdf[['ecdf']](sample.subset.mat)
    names(ecdf.sample.mat) <- colnames(sample.subset.mat)
    return(ecdf.sample.mat)
    })
  
  
  # Append the list containing all the Percentile values and 95% CI values
  ecdf.list[[goi.hugo]] <- tissue.ecdf.list[order(sapply(tissue.ecdf.list, median))]
  ci95.list[[goi.hugo]] <- lapply(tpm.tissue.subset.ecdf[order(sapply(tissue.ecdf.list, median))],
                                  function(x) x[['ci95']])
}



#Generate histology subtype-type index
sra.ord.idx <- sra.dict[with(sra.dict, order(body.site, histology)), ]
sra.ord.idx <- unique(sra.ord.idx[,c("body.site", "histology")])


# Plotting
for(each.goi in names(ecdf.list)){
  print(each.goi)
  pdf(file.path(outdir, paste(each.goi, ".tum2gtexCdf.pdf", sep="")))
  require(scales)
  library(RColorBrewer)
  split.screen(matrix(c(0, 0.8, 0, 0.85,  
                        0.8, 1, 0, 0.85,
                        0, 0.8, 0.85, 1,
                        0.8, 1, 0, 0.85), byrow=TRUE, ncol=4))
  #Tissue colours
  all_col <- unlist(mapply(brewer.pal, c(8, 8, 9, 9), c('Accent', 'Dark2', 'Pastel1', 'Set1')))
  tissue.col <- all_col[c(1:length(all.tissues$V1))]
  names(tissue.col) <- all.tissues$V1
  col.order <- match(names(ecdf.list[[each.goi]]), sra.ord.idx$body.site)
  tissue.col.ord <- tissue.col[sra.ord.idx[col.order, 'histology']]
  
  # Add boxplot
  screen(1)
  par(mar=c(11, 4.1, 0, 0))
  plotTissueBoxplots(ecdf.list[[each.goi]], min.val=0, max.val=1, 
                     tissue.col.ord,  col.n='Percentile')  
  
  # Adds Tissue specific labels
  screen(2)
  par(mar=c(11, 0, 0, 2.1))
  plotTissueLegend(tissue.col)
  
  # Add 95% CI boxes of TPM for each tissue
  screen(3)
  par(mar=c(0, 4.1, 2.1, 0))
  plot95Ci(ci95.list[[each.goi]], min.val=0, max.val=150, 
           tissue.col=tissue.col.ord, goi=each.goi)
  
  close.screen(all.screens=TRUE)
  dev.off()
  
  
}
