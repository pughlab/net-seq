library(reshape2)
library(ggplot2)

# Load in GTEX Data
load("/mnt/work1/users/pughlab/external_data/GTEx/FPKMmatrix/tpm/rdata/Pancreas.annotation.Rdata")
load("/mnt/work1/users/pughlab/external_data/GTEx/FPKMmatrix/tpm/rdata/Pancreas.tpm.Rdata")
gtex.pancreas.mat <- tpm.mat
rm.idx <- which(apply(gtex.pancreas.mat, 2, function(x) all(is.na(x))))
gtex.pancreas.mat <- gtex.pancreas.mat[,-rm.idx]
gtex.pancreas.mat <- normalize.quantiles(gtex.pancreas.mat)
colnames(gtex.pancreas.mat) <- colnames(tpm.mat)[-rm.idx]
rownames(gtex.pancreas.mat) <- rownames(tpm.mat)

# Load in NETseq Data
load("/mnt/work1/users/pughlab/projects/NET-SEQ/loh_signature/wes_validation-cohort/data/loh-het_regions/NETseq_loh-het-genes.Rdata")
load("/mnt/work1/users/pughlab/projects/NET-SEQ/rna_seq/Cufflinks/star2.4.2a/data_matrix/tpm/rdata/NETSEQ_on-treatment.annotation.Rdata")
load("/mnt/work1/users/pughlab/projects/NET-SEQ/rna_seq/Cufflinks/star2.4.2a/data_matrix/tpm/rdata/NETSEQ_pre-treatment.tpm.Rdata")
pre.tpm.mat <- tpm.mat
load("/mnt/work1/users/pughlab/projects/NET-SEQ/rna_seq/Cufflinks/star2.4.2a/data_matrix/tpm/rdata/NETSEQ_on-treatment.tpm.Rdata")
tpm.mat <- cbind(tpm.mat, pre.tpm.mat)

global.anno.df$ensembl <- gsub("\\..*", "", global.anno.df$gene_id)
het.anno.idx <- match(het.genes$ensembl, global.anno.df$ensembl)
loh.anno.idx <- match(loh.genes$ensembl, global.anno.df$ensembl)

pnet.samples <- paste('NET', c('001', '003', '008', '009'), 'T_4', sep='-')
ginet.samples <- paste('NET', c('002', '006', '007'), 'T_4', sep='-')

###############
## Functions ##
###############
gg_color_hue <- function(n) {
  hues = seq(15, 375, length = n + 1)
  hcl(h = hues, l = 65, c = 100)[1:n]
}

avgExpr <- function(tpm.mat, row.idx, col.idx, group.type){
  avg.mat <- matrix(apply(tpm.mat[row.idx, col.idx], 1, mean), ncol=1)
  colnames(avg.mat) <- group.type
  rownames(avg.mat) <- rownames(tpm.mat)[row.idx]
  return(avg.mat)
}

meltExprMat <- function(tpm.mat, row.idx, col.idx=NA, 
                        net.type, loh.type, group.type){
  if(all(is.na(col.idx))) col.idx <- c(1:ncol(tpm.mat))
  if(group.type != 'Sample') {
    tpm.mat <- avgExpr(tpm.mat, row.idx, col.idx, group.type)
    melt.tpm <- melt(tpm.mat)
  } else {
    melt.tpm <- melt(tpm.mat[row.idx, col.idx])
  }
  melt.tpm$net.type <- net.type
  melt.tpm$loh.stat <- loh.type
  melt.tpm$group.stat <- group.type
  colnames(melt.tpm) <- c("gene", "s.id", "value", "net.type", "loh.stat", "group.stat")
  
  melt.tpm$value <- log2(melt.tpm$value + 1)
  return(melt.tpm)
}

getQuant <- function(x, low.perc=0.25, hi.perc=0.75){
  med.val <- median(x$value)
  low.percentile <- quantile(x$value, low.perc)
  hi.percentile <- quantile(x$value, hi.perc)
  return(data.frame("lo"=round(low.percentile,2),
                    "med"=round(med.val,2),
                    "hi"=round(hi.percentile,2)))
}


##########
## Main ##
##########
# Assemble the melted data structure for plotting
melt.list <- list()
melt.list[['pnet-het']] <- meltExprMat(tpm.mat, het.anno.idx, pnet.samples, "PNET", "Het", "Sample")
melt.list[['ginet-het']] <- meltExprMat(tpm.mat, het.anno.idx, ginet.samples, "GINET", "Het", "Sample")
melt.list[['pnet-loh']] <- meltExprMat(tpm.mat, loh.anno.idx, pnet.samples, "PNET", "LOH", "Sample")
melt.list[['ginet-loh']] <- meltExprMat(tpm.mat, loh.anno.idx, ginet.samples, "GINET", "LOH", "Sample")

melt.list[['group-pnet-loh']] <- meltExprMat(tpm.mat, loh.anno.idx, pnet.samples, "PNET", "LOH", "NETseq-PNET")
melt.list[['group-ginet-loh']] <- meltExprMat(tpm.mat, loh.anno.idx, ginet.samples, "GINET", "LOH", "NETseq-GINET")
melt.list[['group-pnet-het']] <- meltExprMat(tpm.mat, het.anno.idx, pnet.samples, "PNET", "Het", "NETseq-PNET")
melt.list[['group-ginet-het']] <- meltExprMat(tpm.mat, het.anno.idx, ginet.samples, "GINET", "Het", "NETseq-GINET")

melt.list[['gtex-panc-loh']] <- meltExprMat(gtex.pancreas.mat, loh.anno.idx, NA, "Pancreas", "LOH", "GTEx")
melt.list[['gtex-panc-het']] <- meltExprMat(gtex.pancreas.mat, het.anno.idx, NA, "Pancreas", "Het", "GTEx")

melt.df <- do.call("rbind", melt.list)
melt.df$s.id <- as.character(melt.df$s.id)
melt.df$s.id <- factor(melt.df$s.id,
                       levels=c(pnet.samples, "NETseq-PNET", "GTEx", ginet.samples, "NETseq-GINET"))


lapply(melt.list, function(xmelt){
  x.sid <- split(xmelt, xmelt$s.id)
  sapply(x.sid, getQuant)
})

# Generate violin/boxplots
pdf("vioplot.pdf", width=15, height = 6)
# TPM for each sample, grouped by LOH/Het status and cohorts
for(each.lohstat in unique(melt.df$loh.stat)){
  melt.loh.df <- melt.df[which(melt.df$loh.stat == each.lohstat), ]
  print(ggplot(melt.loh.df, aes(x = s.id, y = value, fill = net.type)) + 
    geom_violin() + geom_boxplot(width = 0.2) +
    facet_grid(~ group.stat) + 
    theme_bw() + 
    theme(panel.border = element_blank(), panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"),
          axis.text.x = element_text(angle = 90, hjust = 1)))
}


# TPM for each het/loh region, grouped by sample ids
ggplot(melt.df, aes(x = loh.stat, y = value, fill = net.type)) + 
  geom_violin() + geom_boxplot(width = 0.2) +
  facet_grid( ~ s.id) + 
  theme_bw() + 
  theme(panel.border = element_blank(), panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"),
        axis.text.x = element_text(angle = 90, hjust = 1))

dev.off()
