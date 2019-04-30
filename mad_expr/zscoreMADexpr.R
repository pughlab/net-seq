### Z-scores the extended NET data based on MAD status...
### and maybe adds a gistic style summary plot?
library(AnnotationDbi)
library(org.Hs.eg.db)
library(scales)
library(preprocessCore)
plotpdf <- FALSE

#### Load in Data #### 
## Load in MAD samples
mad.dir <- '/mnt/work1/users/pughlab/projects/NET-SEQ/rna_seq_external/MAE'
MAD.samples <- file.path(mad.dir, "samples", "samples_MAD-netseq.txt")
noLOH.samples <- file.path(mad.dir, "samples", "samples_nonLOH.txt")

MAD <- read.table(MAD.samples, header=FALSE,
                  stringsAsFactors = FALSE, check.names = FALSE)
MAD <- as.character(unlist(MAD))
mad.neg <- read.table(noLOH.samples, header=FALSE,
                   stringsAsFactors = FALSE, check.names = FALSE)
mad.neg <- as.character(unlist(mad.neg))
mad.pos <- gsub("(.vcf.annotated)?.maf", "", MAD[!MAD %in% mad.neg])
mad.neg <- gsub("(.vcf.annotated)?.maf", "", mad.neg)

## Load in MAD samples
script.dir <- '/mnt/work1/users/home2/quever/git/loh_cn_visualize/data'
load(file.path(script.dir, "NETseq_loh-het-genes.Rdata"))  # het.henes, loh.genes

## Load in TPM
pdir <- '/mnt/work1/users/pughlab/projects/NET-SEQ/rna_seq_external/clustering'
# Loads in all TPM matrices and records their IDs/class
loaded.data <- lapply(list.files(file.path(pdir, "input"), 
                                 pattern = "tpm.*Rdata"), function(tpm.id){
                                   setwd(pdir)
                                   uid <- gsub(".tpm.*", "", tpm.id)
                                   print(paste0("Loading in ", uid, "..."))
                                   anno.id <- paste0(uid, ".annotation.Rdata")
                                   load(file.path("input", anno.id))
                                   global.anno.df <<- global.anno.df
                                   load(file.path("input", tpm.id))
                                   
                                   tpm.mat <- round(tpm.mat, 2)
                                   id.class <- data.frame(class=rep(uid, ncol(tpm.mat)),
                                                          id=colnames(tpm.mat))
                                   
                                   list("tpm"=tpm.mat, "id"=id.class)
                                 })


# Clean up the TPM matrix:
g.tpm.mat <- do.call("cbind", lapply(loaded.data, function(x) x[['tpm']]))
g.ids <- do.call("rbind", lapply(loaded.data, function(x) x[['id']]))

g.tpm.mat[is.na(g.tpm.mat)] <- 0
t.tpm.mat <- t(g.tpm.mat)
zeroes <- which(apply(t.tpm.mat, 2, quantile, 0.95) == 0)
t.tpm.mat <- t.tpm.mat[,-zeroes]
global.anno.df <- global.anno.df[-zeroes,]

# Add Chr and pos columns for gene annotation
global.anno.df$chr <- factor(with(global.anno.df, gsub(":.*", "", locus)),
                              levels=c(1:22, "X", "Y", "MT"))
global.anno.df$spos <- as.numeric(with(global.anno.df, 
                                       gsub(".*:([0-9]+)-.*", "\\1", locus)))
global.anno.df$epos <- as.numeric(with(global.anno.df, 
                                       gsub("^.*-", "", locus)))

# Reorder based on chromosome and position order
ord <- with(global.anno.df, order(chr, spos))
t.tpm.mat <- t.tpm.mat[,ord]
global.anno.df <- global.anno.df[ord,]

anno.split <- split(global.anno.df, f=global.anno.df$chr)
cum.x <- sapply(anno.split, function(x) c(min(x$spos), diff(x$spos)))
cum.x <- cumsum(as.numeric(unlist(cum.x[c(1:22, "X", "Y", "MT")])))
length(cum.x) <- nrow(global.anno.df)
global.anno.df$cumx <- cum.x

rle.x <- rle(as.character(global.anno.df$chr))
chg <- c(1, cumsum(rle.x$length))
chr.xpos <- global.anno.df$cumx[chg]




#### Analyze the TPM matrices #### 
# madneg_gtex-pancreas, madpos_gtex-pancreas, madpos_madneg
# madneg, gtex-pancreas, madpos
f1='madpos' 
f2='madneg'
file=paste0(f1, "_", f2)
action="load"

if(action=='save'){
  ## Subset for MAD+ and MAD- samples
  if(f1=='madneg'){
    mad.pos.tpm <- t.tpm.mat[which(rownames(t.tpm.mat)  %in% mad.neg), ]
  } else if(f1=='madpos'){
    mad.pos.tpm <- t.tpm.mat[which(rownames(t.tpm.mat)  %in% mad.pos), ]
  } else {
    NULL
  }
  f1.ids <- rownames(mad.pos.tpm)
  
  if(f2=='madneg'){
    mad.neg.tpm <- t.tpm.mat[which(rownames(t.tpm.mat)  %in% mad.neg), ]
  } else if(f2=='gtex-pancreas'){
    mad.neg.tpm <- t.tpm.mat[which(rownames(t.tpm.mat) %in%
                                     as.character(split(g.ids, f=g.ids$class)[['Pancreas']]$id)),]
  } else {
    NULL
  }
  f2.ids <- rownames(mad.neg.tpm)
  
  # Quantile normalize the new set
  m.tpm <- rbind(mad.pos.tpm, mad.neg.tpm)
  qn.m.tpm <- t(normalize.quantiles(t(m.tpm)))
  madpos.idx <- seq_along(f1.ids)
  mad.pos.tpm <- as.data.frame(qn.m.tpm[madpos.idx,])
  mad.neg.tpm <- as.data.frame(qn.m.tpm[-madpos.idx,])
  colnames(mad.pos.tpm)  <- colnames(mad.neg.tpm) <- colnames(t.tpm.mat)
  
  mad.pos.tpm <- log10(mad.pos.tpm + 1)
  mad.neg.tpm <- log10(mad.neg.tpm + 1)
  
  ## Calculate z-scores of MAD+ samples
  z.pos.mat <- sapply(colnames(mad.pos.tpm), function(gene){
    sapply(mad.pos.tpm[,gene], function(i){
      ref <- mad.neg.tpm[,gene]
      (i - mean(ref, na.rm=TRUE)) / sd(ref, na.rm=TRUE)
    })
  })
  save(mad.pos.tpm, mad.neg.tpm, z.pos.mat, file=paste0(file, ".rda"))
} else if (action == 'load'){
  load(paste0(file, ".rda"))  # mad.pos.tpm, mad.neg.tpm, z.pos.mat
}

mean.z <- apply(z.pos.mat, 2, median)
mean.z[is.infinite(mean.z)] <- quantile(mean.z, 0.99, na.rm=TRUE)

#### Visualize z-score matrices ####
expr <- data.frame("z"=mean.z,
                   "idx"=global.anno.df$cumx)
# Fit the Loess Model
loessMod10 <- loess(z ~ idx, data=expr, span=0.05) # 10% smoothing span
smoothed10 <- predict(loessMod10, expr$idx)

min.x <- -3.5
max.x <- 3.5

if(plotpdf){
  pdf(file.path(pdir, "output", paste0(file, ".pdf")), height = 4)
  reverse.idx <- FALSE
  sample.size <- 0.3
  max.x <- quantile(expr$z, 0.99, na.rm=TRUE)
  min.x <- quantile(expr$z, 0.01, na.rm=TRUE)
  
  plot.idx <- sample(1:nrow(expr), size = sample.size * nrow(expr), replace = FALSE)
  expr.pl <- expr[sort(plot.idx),]
  with(expr.pl, plot(x=if(reverse.idx) -1*idx else idx, (z), 
                  pch=16, col=alpha("grey", 0.3),
                  ylim=c(min.x, max.x), las=1,
                  main="MAD+ Expr", xaxt='n', xlab='', 
                  ylab='Expression (z-score)', 
                  type='p'))
  chr.start <- chr.xpos[-length(chr.xpos)]
  chr.end <- chr.xpos[-1][c(TRUE,FALSE)]
  rect(ytop=-10, ybottom = 10, 
       xleft=if(reverse.idx) -1*chr.start[c(TRUE,FALSE)] else chr.start[c(TRUE,FALSE)], 
       xright = if(reverse.idx) -1*chr.end else chr.end, 
       border=FALSE, col=alpha("blue", 0.1))
  lines(if(reverse.idx)-1*(expr$idx) else expr$idx, smoothed10, 
        lwd=1.5, col="black")
  abline(h=0, lty=2, col="black")
  text(x=if(reverse.idx) (-1*(chr.start + 10)) else (chr.start + 10),
       y=rep(c(min.x, min.x+0.01), length(chr.xpos)),
       labels=gsub("^chr", "", rle.x$values),
       adj=0, cex=0.5)
  dev.off()
}

