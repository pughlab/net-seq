library(plyr)
pdir <- '/mnt/work1/users/pughlab/projects/NET-SEQ/rna_seq_external/MAE'
script.dir <- '/mnt/work1/users/home2/quever/git/loh_cn_visualize/data'
MAD.samples <- file.path(pdir, "samples", "samples_MAD-netseq.txt")
noLOH.samples <- file.path(pdir, "samples", "samples_nonLOH.txt")
load(file.path(script.dir, "NETseq_loh-het-genes.Rdata"))

min.snps <- 1
min.depth <- 1
MAD <- read.table(MAD.samples, header=FALSE,
                  stringsAsFactors = FALSE, check.names = FALSE)
MAD <- as.character(unlist(MAD))
ginets <- read.table(noLOH.samples, header=FALSE,
                     stringsAsFactors = FALSE, check.names = FALSE)
ginets <- as.character(unlist(ginets))
pnets <- MAD[!MAD %in% ginets]
setwd(pdir)

all.het.mafs.avg.af <- lapply(list.files(file.path(pdir, "maf")), function(each.file){
  each.maf <- read.table(file.path(pdir, "maf", each.file), header=TRUE, 
                         stringsAsFactors = FALSE, check.names = FALSE,
                         comment.char = "#", sep="\t")
  het.mafs.avg.af <- tryCatch({
    het.maf <- each.maf[which(each.maf$Hugo_Symbol %in% het.genes$genes),]
    het.maf <- het.maf[order(het.maf$Hugo_Symbol),]
    het.maf.per.gene <- split(het.maf, f=het.maf$Hugo_Symbol)
    {
      snps.per.gene <- sapply(het.maf.per.gene, nrow)
      tbl.snps <- t(table(snps.per.gene))
      tbl.snps.per.gene <- t(snps.per.gene)
      
      
      avg.cov.snp.per.gene <- sapply(het.maf.per.gene, function(x) mean(x$read_depth))
      tbl.cov <- t(floor(avg.cov.snp.per.gene))
    }
    
    
    snp.filt <- which(sapply(het.maf.per.gene, nrow) > min.snps)
    cov.filt <- which(sapply(het.maf.per.gene, function(x) all(x$read_depth > min.depth)))
    snp.cov.filt <- intersect(snp.filt, cov.filt)
    
    het.maf.per.gene <- het.maf.per.gene[snp.cov.filt]
    het.mafs.avg.af <- sapply(het.maf.per.gene, function(x) mean(x$allele_frequency))
    het.mafs.avg.af <- t(het.mafs.avg.af)
    list("hets"=het.mafs.avg.af,
         "snps"=tbl.snps,
         "snps.per.gene"=tbl.snps.per.gene,
         "cov"=tbl.cov)
  }, error=function(e) { list("hets"=NA,
                              "snps"=NA,
                              "snps.per.gene"=NA,
                              "cov"=NA) })
  
  het.mafs.avg.af
})

keep.idx <- sapply(all.het.mafs.avg.af, function(x) length(x[['hets']]) > 1)
keep.ids <- list.files(file.path(pdir, "maf"))[which(keep.idx)]
all.het.mafs.avg.af <- all.het.mafs.avg.af[which(keep.idx)]
names(all.het.mafs.avg.af) <- keep.ids

snps.mat <- rbind.fill(lapply(all.het.mafs.avg.af, function(x) as.data.frame.matrix(x[['snps']])))
rownames(snps.mat) <- keep.ids
snps.per.gene.mat <- rbind.fill(lapply(all.het.mafs.avg.af, function(x) as.data.frame.matrix(x[['snps.per.gene']])))
rownames(snps.per.gene.mat) <- keep.ids
cov.mat <- rbind.fill(lapply(all.het.mafs.avg.af, function(x) as.data.frame.matrix(x[['cov']])))
rownames(cov.mat) <- keep.ids
all.het.mat <- rbind.fill(lapply(all.het.mafs.avg.af, function(x) as.data.frame.matrix(x[['hets']])))
rownames(all.het.mat) <- keep.ids




## Working on a per gene basis

pnet.het.mat <- all.het.mat[which(rownames(all.het.mat) %in% pnets),]
ginet.het.mat <- all.het.mat[which(rownames(all.het.mat) %in% ginets),]

na.idx <- apply(pnet.het.mat, 2, function(x) sum(is.na(x)) <= 0)
hom.na.idx <- apply(pnet.het.mat, 2, function(x){
  (sum(is.na(x)) <= 10) && (quantile(na.omit(x),0.05) > 0.8)
} )
pnet.het.mat[,which(na.idx)]
round(pnet.het.mat[,which(hom.na.idx), drop=FALSE],2)



pnet.snps.per.gene <- snps.per.gene.mat[,which(colnames(snps.per.gene.mat) %in% colnames(pnet.het.mat))]
pnet.cov.per.gene <- cov.mat[,which(colnames(cov.mat) %in% colnames(pnet.het.mat))]
## Order all matrices the same way:
match.order.1 <- order(colnames(pnet.snps.per.gene)) == order(colnames(pnet.cov.per.gene))
pnet.het.mat <- pnet.het.mat[,colnames(pnet.snps.per.gene)]
match.order.2 <- order(colnames(pnet.het.mat)) == order(colnames(pnet.cov.per.gene))
if(all(match.order.1) && all(match.order.2)){
  # Take the average weighted by the number of samples with values; order by that
  gene.af <- apply(pnet.het.mat, 2, function(x) mean(x, na.rm=TRUE) * length(which(!is.na(x))))
  gene.ord <- order(gene.af, decreasing = TRUE)
  
  pnet.het.mat <- pnet.het.mat[,gene.ord]
  pnet.snps.per.gene <- pnet.snps.per.gene[,gene.ord]
  pnet.cov.per.gene <- pnet.cov.per.gene[,gene.ord]
  save(pnet.het.mat, pnet.snps.per.gene, pnet.cov.per.gene, 
       file=file.path(pdir, "output", "pnet_het_snps.Rdata"))
}


plotMat <- function(mat.x, xlab, ylab, scr.mar, scr.idx, ylim, max.x=100){
  screen(scr.idx)
  par(mar=scr.mar)
  mat.x[is.na(mat.x)] <- 0
  if(ncol(mat.x) > max.x) mat.x <- mat.x[,1:max.x]
  plot(x=sort(rep(seq(1, ncol(mat.x)), nrow(mat.x))), ylim=c(0, ylim),
       y=unlist(mat.x), ylab=ylab, xlab=xlab, xaxt='n', las=2)
}
non.na.idx <- apply(pnet.het.mat, 2, function(x) length(which(is.na(x))) <= 10)

pdf(file.path(pdir, "output", "pnet-het-retained.pdf"), width=50)

split.screen(c(3, 1))
plotMat(pnet.het.mat[,non.na.idx], "", ylab="Alt Allele Fraction", c(1, 4.1, 4.1, 2.1), 1, 1)
plotMat(pnet.snps.per.gene[,non.na.idx], "", ylab="Snps per gene", c(0.5, 4.1, 0.5, 2.1), 2, 50)
plotMat(pnet.cov.per.gene[,non.na.idx], "", ylab="Cov per gene", c(5.1, 4.1, 0.5, 2.1), 3, 100)
axis(side = 1, at=seq(1:length(which(non.na.idx))),
     labels = colnames(pnet.het.mat)[non.na.idx], las=2, cex.axis=0.5)
close.screen(all.screens=TRUE)

dev.off()
