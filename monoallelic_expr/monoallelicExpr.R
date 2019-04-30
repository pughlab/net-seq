library(weights)
library(plyr)
library(scales)
library(survcomp)
library(VennDiagram)
library(org.Hs.eg.db)
library(TxDb.Hsapiens.UCSC.hg19.knownGene)
library(reshape2)

pdir <- '/mnt/work1/users/pughlab/projects/NET-SEQ/rna_seq_external/MAE'
script.dir <- '/mnt/work1/users/home2/quever/git/loh_cn_visualize/data'
MAD.samples <- file.path(pdir, "samples", "samples_MAD-netseq.txt")
noLOH.samples <- file.path(pdir, "samples", "samples_nonLOH.txt")
load(file.path(script.dir, "NETseq_loh-het-genes.Rdata"))

#### PARAMETERS ####
max.frac.NA.samples <- 0.25  # max fraction of samples to have no read coverage
ordering <- TRUE # boolean whether to order genes based on AF
AF.threshold <- 0.8 # Allelic imbalance threshold where 95% of data must be above
min.snps <- 2
q.cutoff <- 0.2 # Samples must have a q value less than this to be considere "Significant"
min.depth <- 10  # 1 - pbinom(3,4,0.44) is the min to get 0.05 pval
min.samples <- 0.50 # Min fraction of samples to needed to evaluate MAE
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
                                           colnames(pval.mat),
                                           colnames(estq.mat),
                                           colnames(estp.mat)))
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
  out.estp.mat <- estp.mat[r.idx,gene.intersect]
  out.estq.mat <- estq.mat[r.idx,gene.intersect]
  
  
  list("Snps"=out.snps.per.gene,
       "Cov"=out.cov.per.gene,
       "AF"=out.het.mat,
       "p"=out.pval.mat,
       'estp'=out.estp.mat,
       'estq'=out.estq.mat)
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

# Visualization: Plots points on the boxplot and colours according to filters
overlapPoints <- function(melt.af, q, chr.ord,
                          gene.idx, jitter=0.3,
                          print.chr=TRUE, print.stats=FALSE,
                          plot.stats=TRUE,
                          sig.col='red', nonsig.col="gray39", 
                          nonpow.col='black'){
  
  melt.af$xidx <- rep(rev(gene.idx), 
                      rle(as.character(melt.af$gene))$lengths)
  melt.af$xidx.jitter <- melt.af$xidx + runif(n = nrow(melt.af), 
                                              min = (-1 * jitter), max = jitter)
  
  melt.q <- melt(q[,chr.ord])
  melt.af$q <- melt.q$value
  
  melt.af <- melt.af[-grep("^chr", melt.af$gene),]
  # Non-sig q samples per gene
  with(melt.af[which(!is.na(melt.af$q)),],
       points(y=xidx.jitter, x=value, col=alpha(nonsig.col, 0.6), pch=16))
  
  # Underpowered and untested samples/gene
  with(melt.af[which(is.na(melt.af$q)),], 
       points(y=xidx.jitter, x=value, col=alpha(nonpow.col, 1), pch=1))
  
  # Significant q genes/samples
  with(melt.af[which(melt.af$q < 0.2),],
       points(y=xidx.jitter, x=value, col=alpha(sig.col, 0.6), pch=16))
  
  spl.af <- split(melt.af, f=as.character(melt.af$gene))
  summ.n <- sapply(spl.af, function(x) {
    sigq <- length(which(x$q < 0.2))
    lowq <- length(which(with(x, q >= 0.2 & !is.na(value))))
    uncov <- length(which(with(x, is.na(q) & !is.na(value))))
    
    c("sigq"=sigq, "lowq"=lowq, "uncov"=uncov)
  })
  if(print.stats){
    summ <- apply(summ.n, 2, function(i) paste0("[", paste(i, collapse=","), "]"))
    
    text(x=0.4, y = sapply(spl.af, function(x) unique(x$xidx)),
         labels=summ, pos = 4, cex=0.8)
  } else if(plot.stats){
    n <- unique(sapply(spl.af, nrow))
    
    prop <- apply(summ.n, 2, function(i) round(i / n, 2)/10)
    prop.xmin <- as.numeric(apply(prop, 2, function(i) c(0, cumsum(i)[-3])))
    prop.xmax <- as.numeric(apply(prop, 2, function(i) cumsum(i)))
    
    y.idx <- matrix(rep(sapply(spl.af, function(x) unique(x$xidx)),3), ncol=3)
    y.idx <- as.integer(t(y.idx))
    
    cols <- rep(c(alpha(sig.col, 0.6), 
                  alpha(nonsig.col, 0.6), 
                  alpha(nonpow.col, 0)), length(spl.af))
    
    min.x <- 0.35
    rect(xleft = min.x + prop.xmin, ybottom = y.idx - 0.5, 
         xright = min.x + prop.xmax, ytop = y.idx + 0.5,
         col=cols, border=TRUE)
  }
  
  if(print.chr){
    spl.af <- split(melt.af, f=melt.af$chr)
    chr.pos <- sapply(spl.af, function(x) mean(x$xidx))
    
    axis(side = 2, at = c(max(melt.af$xidx) + 1, chr.pos), 
         labels = c("Chr", gsub("chr", "", names(chr.pos), ignore.case = TRUE)), 
         cex=0.7, line=3, tick = FALSE, las=1)
  }
  
  melt.af
}

# Visualization: Melts dataframes of AF and appends relevant data
meltAf <- function(af, g1, chr.ord, do.all=FALSE){
  # Create chrX_Gene IDs, sort by Chr
  colnames(af) <- paste(as.character(seqnames(g1)), 
                        colnames(af), sep="_")
  
  af.ord <- af[,chr.ord]
  
  # Melt and split Chr/Gene into columns
  melt.af <- melt(af.ord)
  chr.gene <- sapply(as.character(melt.af$variable), function(x) strsplit(x, split="_"))
  melt.af <- cbind(melt.af, do.call("rbind", chr.gene))
  colnames(melt.af) <- c("variable", "value", "chr", "gene")
  if(do.all) melt.af$gene <- melt.af$chr
  
  # Reverse order of gene levels for boxplot plotting
  g <- as.character(melt.af$gene)
  melt.af$gene <- factor(g, levels=rev(unique(g)))
  melt.af
}

# Accessor for the PNET/GINET data structure to get AF/q-values/Cov/SNPs
getSampleInfo <- function(x, samples=pnet, type='AF', 
                          summarize=TRUE, ord=NA){
  xy <- lapply(x, function(y){
    if(summarize) {
      dat <- summary(c(samples[[type]][,y], NA)) 
    } else {
      dat <- samples[[type]][,y]
    }
    as.data.frame(as.matrix(dat))
  })
  xy <- do.call("cbind",xy)
  colnames(xy) <- x
  
  if(all(is.na(ord))){
    if(summarize) {
      ord <- order(xy['Mean',])
    } else {
      ord <- order(apply(xy, 2, mean, na.rm=TRUE), decreasing = TRUE)
    }
  } else {
    ord <- match(ord, colnames(xy))
  }
  xy <- xy[,ord]
  xy
}

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
    
    # Removes intronic and UTR genes with low coverage due to isoform mapping issues
    het.maf.per.gene.filt <- lapply(het.maf.per.gene, function(gene){
      low.depth <- with(gene, read_depth < min.depth)
      premrna.snp <- gene$Variant_Classification %in% c("5UTR", "3UTR", "Intron")
      rm.rows <- which(premrna.snp & low.depth)
      if(length(rm.rows) > 0) gene[-rm.rows,] else gene
    })
    het.maf.per.gene <- het.maf.per.gene.filt
    
    # Calculates descriptive stats: number of snps per gene
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
    het.mafs.avg.af <- t(het.mafs.avg.af)  # Uses oncotoators 0.5 and 1.0 annotation
    
    # Calculates a depth based average allelic fraction per gene
    het.mafs.avg.actual.af <- sapply(het.maf.per.gene, function(x){
      cnt <- lapply(x$allelic_depth, function(i) as.numeric(strsplit(i, ",")[[1]]))
      cnt <- as.numeric(sapply(cnt, function(i) i[which.max(i)]))
      af <- round(cnt/x$read_depth, 2)
      mean(af)
    })
    het.mafs.avg.actual.af <- t(het.mafs.avg.actual.af) 
    
    # Calculates continuous AF and binomial p values
    y <- lapply(het.maf.per.gene, function(x, filter=NA, verbose=FALSE){
      # Get alt/ref allele depth
      cnt <- lapply(x$allelic_depth, function(i) as.numeric(strsplit(i, ",")[[1]]))
      cnt <- as.numeric(sapply(cnt, function(i) i[which.max(i)]))
      af <- round(cnt/x$read_depth, 2)
      if(any(is.na(cnt))){
        cnt[is.na(cnt)] <- 0
        af[is.na(af)] <- 1
      } 

      res <- data.frame("class"=as.character(x$Variant_Classification),
                        "alt"=cnt, "depth"=as.numeric(x$read_depth),
                        "af"=af, stringsAsFactors = FALSE)
      res$P <- round(with(res, (depth-alt) / depth), 3)
      if(any(is.na(res$P))) res <- res[-which(is.na(res$P)),]
      if(any(res$P == 0)) res$P[res$P == 0] <- 10^-(res$depth[res$P == 0]/10)
      if(any(res$P <= 0.001)) res$P[res$P <= 0.001] <- 0.001
      res$RR <- with(res, P/0.5)
      if(any(na.omit(res$RR) > 1)) res$RR[res$RR > 1] <- 1
      res$logRR <- with(res, (-1 * log10(RR)) + 1) 
      
      res
      
    }, filter=NA, verbose=TRUE)
    
    
    null.t <- function(x,y,R=100){
      xy <- round(c(x,y),3)
      x <- (x - mean(x, na.rm=TRUE)) + mean(xy, na.rm=TRUE)
      y <- (y - mean(y, na.rm=TRUE)) + mean(xy, na.rm=TRUE)
      
      sapply(1:R, function(i){
        x0 <- sample(xy, size=length(x), replace = TRUE) 
        y0 <- sample(xy, size=length(y), replace = TRUE) 
        
        round(t.test(x0, y0)$statistic,3)
      })
    }
    
    wtd.t <- function(x,w,y, R=100, min.x=8){
      if(length(x) < min.x) return(NA)
      x <- rep(x, w)
      
      sapply(1:R, function(i){
        x1 <- sample(x, size=length(x), replace = TRUE) 
        y1 <- sample(y, size=length(y), replace = TRUE) 
        
        round(t.test(x1, y1)$statistic,3)
      })
    }
    
    ref.dist <- as.numeric(unlist(sapply(y, function(x) x$af)))
    t0 <- null.t(sample(ref.dist, mean(sapply(y, nrow))), ref.dist, R=1000)
    z <- lapply(y, function(x) wtd.t(x$af, x$logRR, ref.dist))
    est.p <- sapply(z, function(i) sum(t0 > i) / length(t0))
    est.q <- p.adjust(est.p, method = 'fdr', n = length(na.omit(est.p)))
    
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
    
    
    list("hets"=het.mafs.avg.actual.af,
         "rough.hets"=het.mafs.avg.af,
         "snps"=tbl.snps,
         "snps.per.gene"=tbl.snps.per.gene,
         "cov"=tbl.cov,
         "p"=binom.p,
         "est.p"=est.p,
         "est.q"=est.q,
         "pos"=mafs.pos,
         "snp.ks"=y)
  }, error=function(e) { list("hets"=NA,
                              "rough.hets"=NA,
                              "snps"=NA,
                              "snps.per.gene"=NA,
                              "cov"=NA,
                              "p"=NA,
                              "est.p"=NA,
                              "est.q"=NA,
                              "pos"=NA,
                              "snp.ks"=NA) })
  
  het.mafs.avg.af
})

#### Filter Genes ####
# Genes require at least 1 snp covering to be post-filtered
keep.idx <- sapply(all.het.mafs.avg.af, function(x) length(x[['hets']]) > 1)
table(keep.idx)  # 38 keeps, 0 throw-aways
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
# Est.p x Gene matrix
estp.mat <- rbind.fill(lapply(all.het.mafs.avg.af, function(x) as.data.frame(t(x[['est.p']]))))
rownames(estp.mat) <- keep.ids
# Est.q x Gene matrix
estq.mat <- rbind.fill(lapply(all.het.mafs.avg.af, function(x) as.data.frame(t(x[['est.q']]))))
rownames(estq.mat) <- keep.ids


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

#### Summarize significant genes ####
## Overlap KS genes
sig.estq.genes <- apply(pnet[['estq']], 1, function(x) names(which(x < q.cutoff)))
n.cutoff <- nrow(pnet[['estq']]) / 2
sig.estq.genes <- names(which(sort(table(unlist(sig.estq.genes))) > n.cutoff))

g.estq.genes <- apply(ginet[['estq']], 1, function(x) names(which(x < q.cutoff)))
n.cutoff <- nrow(ginet[['estq']]) / 2
g.estq.genes <- names(which(sort(table(unlist(g.estq.genes))) > n.cutoff))



#cat(paste(setdiff(sig.estq.genes, g.estq.genes), collapse="\n"))
all.af.mad <- getSampleInfo(colnames(pnet[['AF']]), 
                            pnet, 'AF', summarize=FALSE)
all.af.wt <- getSampleInfo(colnames(pnet[['AF']]), 
                           ginet, 'AF', summarize=FALSE)

af.mad <- getSampleInfo(setdiff(sig.estq.genes, g.estq.genes), 
                        pnet, 'AF', summarize=FALSE)
q.mad <- getSampleInfo(setdiff(sig.estq.genes, g.estq.genes), 
                       pnet, 'estq', summarize=FALSE,
                       ord=colnames(af.mad))
af.wt <- getSampleInfo(setdiff(sig.estq.genes, g.estq.genes), 
                       ginet, 'AF', summarize=FALSE,
                       ord=colnames(af.mad))
q.wt <- getSampleInfo(setdiff(sig.estq.genes, g.estq.genes), 
                      ginet, 'estq', summarize=FALSE,
                      ord=colnames(af.mad))

#### Visualization ####
input.id <- colnames(af.mad)
entrez.id <- mapIds(org.Hs.eg.db, keys=input.id,
                    column="ENTREZID", keytype="SYMBOL",
                    multiVals="first")
g0 <- genes(TxDb.Hsapiens.UCSC.hg19.knownGene)
g1 <- g0[match(entrez.id, g0$gene_id),]

## Process all the genes
all.input.id <- colnames(all.af.mad)
all.entrez.id <- mapIds(org.Hs.eg.db, keys=all.input.id,
                    column="ENTREZID", keytype="SYMBOL",
                    multiVals="first")
m.idx <- match(all.entrez.id, g0$gene_id)
all.na.idx <- which(is.na(m.idx))
g.all <- g0[m.idx[-all.na.idx],]
all.af.mad <- all.af.mad[,-all.na.idx]
all.af.wt <- all.af.wt[,-all.na.idx]

chr.ord <- order(seqnames(g1))
split.idx <- as.character(seqnames(g1)[chr.ord])

gene.idx <- seq_along(as.character(seqnames(g1)[chr.ord]))
cum.len <- cumsum(seqnames(g1)[chr.ord]@lengths)
for(l in cum.len[-length(cum.len)]){
  gene.idx[(l+1):length(gene.idx)] <- gene.idx[(l+1):length(gene.idx)] + 1
}
gene.idx <- sort(max(gene.idx) - gene.idx) + 1


melt.af.mad <- meltAf(af.mad, g1, chr.ord)
melt.af.wt <- meltAf(af.wt, g1, chr.ord)
melt.all.af.mad <- meltAf(all.af.mad, g.all, 
                          order(seqnames(g.all)), do.all=TRUE)
melt.all.af.wt <- meltAf(all.af.wt, g.all,
                         order(seqnames(g.all)), do.all=TRUE)
melt.af.mad <- combineAf(melt.af.mad, melt.all.af.mad)
melt.af.wt <- combineAf(melt.af.wt, melt.all.af.wt)

gene.idx <- seq_along(levels(melt.af.mad$gene))
for(i in grep("chr", levels(melt.af.mad$gene))[-1]){
  gene.idx[i:length(gene.idx)] <- gene.idx[i:length(gene.idx)] + 1
}

combineAf <- function(af, af.all){
  #af.all <- melt.all.af.mad
  pre.levels <- levels(af$gene)
  
  # Key -> Value dictionary of gene to chromosome
  .dict <- strsplit(x=unique(as.character(af$variable)), split = "_")
  gene.dict <- lapply(.dict, function(i) i[1])
  names(gene.dict) <- sapply(.dict, function(i) i[2])
  
  # Extract relevant chromosomes from all Genes
  keep.idx <- which(af.all$chr %in% sapply(gene.dict, function(i) i ))
  af.all <- af.all[keep.idx,]
  
  # Combine melted structures and re-level
  af.rbind <- rbind(af, af.all)
  
  
  rle.prelvl <- rle(as.character(gene.dict[pre.levels]))
  idx <- cumsum(rle.prelvl$lengths)
  names(idx) <- rle.prelvl$values
  
  splitAt <- function(x, pos) unname(split(x, cumsum(seq_along(x) %in% pos)))
  for(each.idx in seq_along(idx[-1])){
    spl.lvl <- splitAt(pre.levels, rev(idx)[each.idx + 1]+1)
    pre.levels <- c(spl.lvl[[1]], 
                    names(rev(idx)[each.idx]), 
                    spl.lvl[[2]])
  }
  pre.levels <- c(names(idx)[1], pre.levels)
  
  g <- as.character(af.rbind$gene)
  af.rbind$gene <- factor(g, levels=pre.levels)
  af.rbind <- af.rbind[order(af.rbind$gene, decreasing = TRUE),]
  af.rbind
}


# Boxplot per chromosome of AF for q-value selected genes
pdf(file.path(pdir, "output", "MAE_qgenes.pdf"))
split.screen(matrix(c(0, 0.55, 0, 1,
                      0.55, 1.0, 0, 1), byrow = TRUE, ncol=4))
screen(1); par(mar=c(5.1, 5.5, 5.1, 0.3))
boxplot(melt.af.mad$value ~ melt.af.mad$gene, 
        at = gene.idx, ylim=c(0.37, 1.0),
        col = alpha("grey", 0.7), horizontal=TRUE,
        xaxs = FALSE, las=1, outline=FALSE,
        main="MAD+", xlab="AF", cex.axis=0.7, xaxt='n')
axis(side = 1, seq(0.5, 1, by=0.1), labels=seq(0.5, 1, by=0.1), cex.axis=0.7)
axis(side = 1, c(0.35, 0.45), line=-1, tick = FALSE,
     labels=c(0, min(table(melt.af.mad$gene))), cex.axis=0.7)
chr.dividers <- setdiff(c(1:max(gene.idx)), gene.idx)
abline(h = chr.dividers, lty=2, col="grey")
abline(v = 0.45)

# Add points and colour by filtered/used
overlapPoints(melt.af.mad, q.mad, 
              chr.ord, gene.idx, jitter=0.3)

screen(2); par(mar=c(5.1, 0.3, 5.1, 2.1))
boxplot(melt.af.wt$value ~ melt.af.wt$gene, 
        at = gene.idx, ylim=c(0.37, 1.0),
        col = alpha("grey", 0.7), horizontal=TRUE,
        xaxs = FALSE, las=1, outline=FALSE, yaxt='n',
        main="MAD-", xlab="AF", cex.axis=0.7, xaxt='n')
axis(side = 1, seq(0.5, 1, by=0.1), labels=seq(0.5, 1, by=0.1), cex.axis=0.7)
axis(side = 1, c(0.35, 0.45), line=-1, tick = FALSE,
     labels=c(0, min(table(melt.af.wt$gene))), cex.axis=0.7)
chr.dividers <- setdiff(c(1:max(gene.idx)), gene.idx)
abline(h = chr.dividers, lty=2, col="grey")
abline(v = 0.45)

# Add points and colour by filtered/used
overlapPoints(melt.af.wt, q.wt, 
              chr.ord, gene.idx, jitter=0.3,
              print.chr=FALSE)
close.screen(all.screens = TRUE)
dev.off()

