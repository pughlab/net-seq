library(plyr)
library(scales)
library(survcomp)
pdir <- '/mnt/work1/users/pughlab/projects/NET-SEQ/rna_seq_external/MAE'
script.dir <- '/mnt/work1/users/home2/quever/git/loh_cn_visualize/data'
MAD.samples <- file.path(pdir, "samples", "samples_MAD-netseq.txt")
noLOH.samples <- file.path(pdir, "samples", "samples_nonLOH.txt")
load(file.path(script.dir, "NETseq_loh-het-genes.Rdata"))

#### PARAMETERS ####
max.frac.NA.samples <- 0.25  # max fraction of samples to have no read coverage
ordering <- TRUE # boolean whether to order genes based on AF
AF.threshold <- 0.8 # Allelic imbalance threshold where 95% of data must be above
min.snps <- 1
min.samples <- 0.50 # Min fraction of samples to needed to evaluate MAE
min.depth <- 4  # 1 - pbinom(3,4,0.44) is the min to get 0.05 pval
mad.multiplier=4 # Removes SNPs with coverage > median + 4MAD per gene
exp.hom.frac=0.44 # Estimate based on  https://doi.org/10.1371/journal.pgen.1000160
use.af=TRUE  # Assigns HOM/HET based on alleic fractions
af.thresh=0.75  # HOM are SNPs with a AF > this threshold
chr.size <- (249*10^6) # For plotting, size of Chr1

#### FUNCTIONS ####
plotMat <- function(mat.x, xlab, ylab, scr.mar, scr.idx, ylim, max.x=100){
  screen(scr.idx)
  par(mar=scr.mar)
  mat.x[is.na(mat.x)] <- 0
  if(ncol(mat.x) > max.x) mat.x <- mat.x[,1:max.x]
  plot(x=sort(rep(seq(1, ncol(mat.x)), nrow(mat.x))), ylim=c(0, ylim),
       y=unlist(mat.x), ylab=ylab, xlab=xlab, xaxt='n', las=2)
}

# Using oncotatot format, takes in a dataframe of a single gene and
# calculates the binomial probability of allelic imbalance
# for one sample
getBinomialP <- function(gene.maf, min.depth=4, 
                         mad.multiplier=4, exp.hom.frac=0.55,
                         use.af=TRUE, af.thresh=0.75){
  # Remove SNPs that are under-powered or artifactual
  depth <- gene.maf$depth_across_samples
  max.depth = median(depth) + mad.multiplier*mad(depth)
  depth.filter <- which((depth > min.depth) & 
                          (depth < max.depth))  
  
  # Set Hom/Het status based on allelic fraction
  if(use.af){
    af <- do.call("rbind", strsplit(gene.maf$allelic_depth, split=","))
    class(af) <- 'numeric'
    af <- round(apply(af / rowSums(af), 1, max),2)
    hom.idx = af > af.thresh
    gene.maf$allele_frequency[hom.idx] <- 1.0
    gene.maf$allele_frequency[!hom.idx] <- 0.5
  }
  
  # binomial test of homozygous snp counts
  if(length(depth.filter) > 0){
    gene.maf <- gene.maf[depth.filter,]
    hom.cnt <- sum(gene.maf$allele_frequency == 1.0)
    tot.snps = nrow(gene.maf)
    p = 1 - pbinom(hom.cnt - 1, tot.snps, exp.hom.frac)
  } else {
    p = NULL
  }
  
  return(p)
}

# Reduces a p-value matrix by using a z-transform
combine.pval <- function(p.mat, resample='resample', resampling.quantile=0.99){
  sampled.pval <- apply(p.mat, 2, function(p){
    p <- as.numeric(p)[!is.na(p)]
    if(resample=='resample'){
      sapply(1:50, function(x){
        pi <- sample(p, size=ceiling(length(p) * 0.75), 
                     replace = FALSE)
        combine.test(pi, method = "z.transform", na.rm = TRUE)
      })
    } else if(resample=='bootstrap'){
      require(boot)
      stop("Bootstrapping is still under development and will not work yet")
      comb <- function(data, indices){
        combine.test(data[indices], method = "z.transform", na.rm = TRUE)
      }
      boot.ci(boot(data=p, statistic=comb, R=1000))
    }
  })
  aggregate.pval <- apply(sampled.pval, 2, quantile, 
                          probs=resampling.quantile, na.rm=TRUE)
  aggregate.pval
}


# Given a matrix of imbalance p-values and the AF per gene, this calculates the
# combined p-value and FDR adjustment
allelicImbalance <- function(p.mat, af.mat, q=0.05){
  aggregate.pval <- combine.pval(p.mat, resample='resample', resampling.quantile=0.99)
  qval <- p.adjust(aggregate.pval, n = length(aggregate.pval),
                   method='fdr')
  
  #From the significant gene list, remove genes where not
  # enough samples had a value
  q.genes <- sort(qval[which(qval <= q)], decreasing=TRUE)
  non.na.cnt <- apply(p.mat[,names(q.genes)], 2, function(x) sum(!is.na(x)))
  powered.genes <- non.na.cnt > (min.samples * nrow(p.mat))
  q.genes <- q.genes[which(powered.genes)]
  # Manual inspection
  q.genes.pval <- round(p.mat[,names(q.genes)], 3)
  q.genes.af <- round(af.mat[,names(q.genes)], 3)
  
  
  list("q"=qval,
       "sig.q"=q.genes,
       "sig.pvals"=q.genes.pval,
       "sig.af"=q.genes.af)
}

# Extracts matrices for SNPs_per_gene, Avg_Cov_per_gene,
# Avg_AF_per_gene and pval_per_gene
mapMatrices <- function(het.mat, r.idx){
  gene.intersect <- Reduce(intersect, list(colnames(all.pos.mat),
                                           colnames(snps.per.gene.mat),
                                           colnames(het.mat),
                                           colnames(cov.mat),
                                           colnames(pval.mat)))
  gene.pos.mat <- all.pos.mat[,gene.intersect]
  genomic.order <- order(as.numeric(gene.pos.mat['chr',]),
                         as.numeric(gene.pos.mat['start.pos',]))
  gene.pos.mat <- gene.pos.mat[,genomic.order]
  gene.intersect <- colnames(gene.pos.mat)
  
  # Subset for rows of interest
  out.snps.per.gene <- snps.per.gene.mat[r.idx,gene.intersect]
  out.cov.per.gene <- cov.mat[r.idx,gene.intersect]
  out.het.mat <- het.mat[r.idx,gene.intersect]
  out.pval.mat <- pval.mat[r.idx,gene.intersect]
  
  list("Snps"=out.snps.per.gene,
       "Cov"=out.cov.per.gene,
       "AF"=out.het.mat,
       "p"=out.pval.mat)
}

# Given  a matrix of AF with columns as genes, it calculates the mean AF
# weighted by the number of samples that have coverage
rankWeightedAF <- function(af.mat, p=1, top.q=0.95){
  # Weight = number of samples with non NA gene expression
  gene.af <- apply(af.mat, 2, function(x){
    mean(x, na.rm=TRUE) * (p * length(which(!is.na(x))))
  })
  order.af <- order(gene.af, decreasing = TRUE)
  top.genes <- sort(gene.af[which(gene.af > quantile(gene.af, top.q, na.rm=TRUE))], 
                    decreasing = TRUE)
  
  list("af"=gene.af, "order"=order.af, "top"=top.genes)
}

# Print a vector for copy/paste into Enrichr
catGenes <- function(x) { cat(paste0(x, "\n")) }

#### Load Data ####
MAD <- read.table(MAD.samples, header=FALSE,
                  stringsAsFactors = FALSE, check.names = FALSE)
MAD <- as.character(unlist(MAD))
ginets <- read.table(noLOH.samples, header=FALSE,
                     stringsAsFactors = FALSE, check.names = FALSE)
ginets <- as.character(unlist(ginets))
pnets <- MAD[!MAD %in% ginets]
setwd(pdir)

#### Parse MAF File ####
all.het.mafs.avg.af <- lapply(list.files(file.path(pdir, "maf")), 
                              function(each.file){
  each.maf <- read.table(file.path(pdir, "maf", each.file), header=TRUE, 
                         stringsAsFactors = FALSE, check.names = FALSE,
                         comment.char = "#", sep="\t")
  het.mafs.avg.af <- tryCatch({
    # Isolates for Genes of Interest, splits based on genes
    het.maf <- each.maf[which(each.maf$Hugo_Symbol %in% het.genes$genes),]
    het.maf <- het.maf[order(het.maf$Hugo_Symbol),]
    het.maf.per.gene <- split(het.maf, f=het.maf$Hugo_Symbol)
    
    # Calculates descriptivestats: number of snps per gene
    snps.per.gene <- sapply(het.maf.per.gene, nrow)
    tbl.snps <- t(table(snps.per.gene))
    tbl.snps.per.gene <- t(snps.per.gene)
    
    # Calculstes the average/median coverage per SNP per gene
    avg.cov.snp.per.gene <- sapply(het.maf.per.gene, function(x) median(x$read_depth))
    tbl.cov <- t(floor(avg.cov.snp.per.gene))
    
    # Flags SNPs with too few SNPs or poor coverage
    snp.filt <- which(sapply(het.maf.per.gene, nrow) > min.snps)
    cov.filt <- which(sapply(het.maf.per.gene, function(x) any(x$read_depth > min.depth)))
    snp.cov.filt <- intersect(snp.filt, cov.filt)
    
    # Calculates a rough average allelic fraction per gene
    het.maf.per.gene <- het.maf.per.gene[snp.cov.filt]
    het.mafs.avg.af <- sapply(het.maf.per.gene, function(x) mean(x$allele_frequency))
    het.mafs.avg.af <- t(het.mafs.avg.af)
    
    # Calculates binomial probability of Homozygous skew in AF
    binom.p <- sapply(het.maf.per.gene, getBinomialP, 
                      min.depth, mad.multiplier, exp.hom.frac, 
                      use.af, af.thresh)
    binom.p <- unlist(binom.p)
  
    ##snp metadata and positions (range)
    mafs.pos <- lapply(het.maf.per.gene, function(x) {
      data.frame(gene=unique(as.character(x$Hugo_Symbol)),
                 chr=unique(as.character(x$Chromosome)),
                 start.pos=min(x$Start_position),
                 end.pos=max(x$End_position), 
                 stringsAsFactors = FALSE)
    })
    mafs.pos <- t(do.call("rbind", mafs.pos))
    
    
    list("hets"=het.mafs.avg.af,
         "snps"=tbl.snps,
         "snps.per.gene"=tbl.snps.per.gene,
         "cov"=tbl.cov,
         "p"=binom.p,
         "pos"=mafs.pos)
  }, error=function(e) { list("hets"=NA,
                              "snps"=NA,
                              "snps.per.gene"=NA,
                              "cov"=NA,
                              "p"=NA,
                              "pos"=NA) })
  
  het.mafs.avg.af
})

#### Filter Genes ####
# Genes require at least 1 snp covering to be post-filtered
keep.idx <- sapply(all.het.mafs.avg.af, function(x) length(x[['hets']]) > 1)
keep.ids <- list.files(file.path(pdir, "maf"))[which(keep.idx)]
all.het.mafs.avg.af <- all.het.mafs.avg.af[which(keep.idx)]
names(all.het.mafs.avg.af) <- keep.ids

#### Assemble (Sample x Gene) matrix ####
# Sample x #Snps/gene matrix: Count of how many genes have 1 snp/gene, etc
snps.mat <- rbind.fill(lapply(all.het.mafs.avg.af, function(x) as.data.frame.matrix(x[['snps']])))
rownames(snps.mat) <- keep.ids
# Sample x Gene matrix: Count of Snps/gene
snps.per.gene.mat <- rbind.fill(lapply(all.het.mafs.avg.af, function(x) as.data.frame.matrix(x[['snps.per.gene']])))
rownames(snps.per.gene.mat) <- keep.ids
# Sample x Gene matrix: Mean coverage across all Snps/gene
cov.mat <- rbind.fill(lapply(all.het.mafs.avg.af, function(x) as.data.frame.matrix(x[['cov']])))
rownames(cov.mat) <- keep.ids
# Sample x Gene matrix: Mean AF across all Snps/gene
all.het.mat <- rbind.fill(lapply(all.het.mafs.avg.af, function(x) as.data.frame.matrix(x[['hets']])))
rownames(all.het.mat) <- keep.ids
# Sample x Gene matrix: binomial pval per gene
pval.mat <- rbind.fill(lapply(all.het.mafs.avg.af, function(x) as.data.frame(t(x[['p']]))))
rownames(pval.mat) <- keep.ids
# Pos x Gene matrix
all.pos.mat <- do.call("cbind", lapply(all.het.mafs.avg.af, function(x) x[['pos']]))
all.pos.mat <- all.pos.mat[,!duplicated(all.pos.mat[1,])]


#### Identify allelic imbalances ####
# Split into meaningful groups, with PNETs being the target 
# samples of interest and ginets are controls
pnet.idx <- which(rownames(all.het.mat) %in% pnets)
ginet.idx <- which(rownames(all.het.mat) %in% ginets)

# Flag Allelicc imbalancebased on high AF
max.num.NA.samples <- ceiling(nrow(all.het.mat) * max.frac.NA.samples)
# hom.na.idx <- apply(pnet.het.mat, 2, function(x){
#   (sum(is.na(x)) <= max.num.NA.samples) && 
#     (quantile(na.omit(x),0.05) > AF.threshold)
# } )
# pnet.het.mat[,which(na.idx)]
# round(pnet.het.mat[,which(hom.na.idx), drop=FALSE],2)



#### Consistency Check ####
# Making sure the 3 matrices, snps/gene, cov/gene, af/gene are all in the same order
pnet <- mapMatrices(all.het.mat, pnet.idx)
ginet <- mapMatrices(all.het.mat, ginet.idx)


## Get significantly imbalanced genes and apply FDR correction
pnet.q <- allelicImbalance(p.mat=pnet[['p']], af.mat=pnet[['AF']], q=0.25)
ginet.q <- allelicImbalance(p.mat=ginet[['p']], af.mat=ginet[['AF']], q=0.25)

## Overlap imbalanced genes with weighted mean AF rank
p.rank.af <- rankWeightedAF(pnet[['AF']], p=1, top=0.95)
g.rank.af <- rankWeightedAF(ginet[['AF']], p=1)

p.sigq <- intersect(names(pnet.q[['sig.q']]), names(p.rank.af[['top']]))
g.sigq <- intersect(names(ginet.q[['sig.q']]), names(g.rank.af[['top']]))

pdf(file.path(pdir, "output", "overlap_pnetMAD.pdf"))
v <- venn.diagram(list(mad.p=p.sigq, mad.n=g.sigq),
                  fill = c("orange", "blue"),
                  alpha = c(0.5, 0.5), cat.cex = 1.5, cex=1.5,
                  filename=NULL)
grid.newpage()
grid.draw(v)
dev.off()

catGenes(setdiff(p.sigq, g.sigq)) ## MAD+ Genes
catGenes(setdiff(g.sigq, p.sigq)) ## MAD- Genes
catGenes(intersect(p.sigq, g.sigq)) ## MAD indepedent genes

#### OPTIONAL/VISUALIZATION: Ordered by most unbalanced AF ####
if(ordering){
  # Take the average weighted by the number of samples with values; order by that
  save(all.het.mafs.avg.af,
       file=file.path(pdir, "output", "pnet_het_snps.Rdata"))
  
  for(ord.type in c('rank', 'sig', 'rank.sig')){
    switch(ord.type, 
           rank=gene.ord <- p.rank.af[['order']],  # All rank-based
           sig=gene.ord <-  match(names(pnet.q[['sig.q']]), colnames(pnet[['AF']])),  # All significant/powered
           rank.sig=gene.ord <-  match(p.sigq, colnames(pnet[['AF']]))  # Rank based + significant/powered
    )
    gene.ord <- head(gene.ord, 50)
    
    pnet.het.mat <- pnet[['AF']][,gene.ord]
    pnet.snps.per.gene <- pnet[['Snps']][,gene.ord]
    pnet.cov.per.gene <- pnet[['Cov']][,gene.ord]
    non.na.idx <- apply(pnet.het.mat, 2, function(x) sum(is.na(x)) <= max.num.NA.samples)

    pdf(file.path(pdir, "output", paste0("pnet-het-retained.", ord.type, ".pdf")), width=15)
    split.screen(c(3, 1))
    plotMat(pnet.het.mat[,non.na.idx], "", ylab="Alt Allele Fraction", c(1, 4.1, 4.1, 2.1), 1, 1)
    plotMat(pnet.snps.per.gene[,non.na.idx], "", ylab="Snps per gene", c(0.5, 4.1, 0.5, 2.1), 2, 50)
    plotMat(pnet.cov.per.gene[,non.na.idx], "", ylab="Cov per gene", c(5.1, 4.1, 0.5, 2.1), 3, 100)
    axis(side = 1, at=seq(1:length(which(non.na.idx))),
         labels = colnames(pnet.het.mat)[non.na.idx], las=2, cex.axis=0.5)
    close.screen(all.screens=TRUE)
    dev.off()
  }
  
}

#### VISUALIZATION: Ordered by genomic coord context ####
non.na.idx <- apply(pnet[['AF']], 2, function(x) sum(is.na(x)) <= max.num.NA.samples)

to.df <- function(x){ data.frame(as.matrix(t(x)), stringsAsFactors=FALSE)  }
pos.mat <- to.df(all.pos.mat)
pos.split <- split(pos.mat, f=pos.mat$chr)

splitNetMatrix <- function(mat, ref, f){
  df <- to.df(mat)
  df <- df[match(ref, rownames(df)),]
  rownames(df) <- ref
  
  split(df[match(ref, rownames(df)),], f)
}
pval.split <- splitNetMatrix(mat=pnet[['p']], ref=rownames(pos.mat), f=pos.mat$chr)
af.split <- splitNetMatrix(mat=pnet[['AF']], ref=rownames(pos.mat), f=pos.mat$chr)
snps.split <- splitNetMatrix(mat=pnet[['Snps']], ref=rownames(pos.mat), f=pos.mat$chr)
cov.split <- splitNetMatrix(mat=pnet[['Cov']], ref=rownames(pos.mat), f=pos.mat$chr)

pdf(file.path(pdir, "output", "pnet-het-retained.genomic.pdf"), height=20, width=15)
chr.ids <- names(af.split)[order(as.numeric(names(af.split)))]
split.screen(c(length(chr.ids), 1))
lapply(chr.ids, function(each.chr){
  pos <- as.numeric(as.character(pos.split[[each.chr]]$start.pos))
  mu <- apply(af.split[[each.chr]], 1, mean, na.rm=TRUE)  # Mean of gene AF
  #mu <- apply(af.split[[each.chr]], 1, median, na.rm=TRUE)  # Median of gene AF
  sd <- apply(af.split[[each.chr]], 1, sd, na.rm=TRUE)
  sd[is.na(sd)] <- 0.5
  n <- apply(pval.split[[each.chr]], 1, function(x) sum(!is.na(x))) 
  n.frac <- n / ncol(pval.split[[each.chr]])
  alpha.frac <- (0.5-sd)*2 * n.frac
  
  screen(grep(paste0("^", each.chr, "$") , chr.ids)); par(mar=c(0.1, 4.1, 0.1, 2.1))
  plot(x=0, y=0, type='n', xaxt='n',
       xlim=c(1, chr.size), ylim=c(0.5,1),  las=1,
       xlab="", ylab=paste0("Chr", each.chr))
  basal.mu <- (0.5 * exp.hom.frac + 0.5)
  abline(h = basal.mu, lty=2, col="grey")
  segments(x0 = pos, y0 = basal.mu, x1 = pos, y1 = mu,
           col=alpha("grey", (alpha.frac)))
  points(x=pos, y=mu, pch=16, cex=(1 + alpha.frac^2),
         col=alpha("gray20", (alpha.frac)))

  psig.mu <- p.sigq[p.sigq %in% names(mu)]
  qval.chr <- pnet.q[['sig.q']][psig.mu]
  q.pos <- match(psig.mu, names(mu))
  if(length(q.pos)>0){
    ord <- order(pos[q.pos])
    points(x=pos[q.pos][ord], y=mu[q.pos][ord], pch=16, cex=(1 + alpha.frac[ord]^2),
           col=alpha("red", (1-qval.chr[ord])))
    text(x=pos[q.pos][ord], y=rep(c(0.48, 0.9), length(q.pos))[1:length(q.pos)], 
         labels = names(mu)[q.pos][ord], pos = 3, cex=0.8, col="red")
  }
})
close.screen(all.screens = TRUE)
dev.off()





