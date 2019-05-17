library(GenomicRanges)
library(PLTK)
library(dplyr)
library(RColorBrewer)
library(scales)
library(metap)
library(mclust)

###################
#### FUNCTIONS ####
#Preprocessing functions
readInWigs <- function(f){
  wig <- read.table(f, header=FALSE, stringsAsFactors = FALSE, 
                    check.names=FALSE, sep="\t")
  colnames(wig) <- c("chr", "start", "end", "peak")
  wig$peak <- round(wig$peak, 1)
  makeGRangesFromDataFrame(wig, keep.extra.columns = TRUE)
}
overlapGrWithROI <- function(r, gr){
  print(r)
  # Isolate for "CEN" regions
  ov <- findOverlaps(gr, cb[which(cb$CEN == r)])
  gr.ov <- gr[queryHits(ov),]
  gr.chrs <- lapply(chrs, function(chr){
    gr.chr <- gr.ov[seqnames(gr.ov)==chr]
    
    # Isolate for groups
    colids <- colnames(elementMetadata(gr.chr))
    col.idx <- sapply(grps, function(grp){grep(grp, colids)})
    
    tryCatch({
      gr.chr$tstat <- apply(elementMetadata(gr.chr), 1, function(i){
        t.test((i[col.idx[,2]]), 
               (i[col.idx[,1]]), 
               alternative = 'less')$statistic
      })
      gr.chr$pval <- apply(elementMetadata(gr.chr), 1, function(i){
        t.test((i[col.idx[,2]]), 
               (i[col.idx[,1]]), 
               alternative = 'less')$p.value
      })
    }, error=function(e){NULL})
    
    gr.chr
  })
  names(gr.chrs) <- chrs
  gr.chrs
}


getCENidx <- function(cb, spread=1){
  cen.idx <- which(cb$gieStain == 'acen')
  pcen.idx <- cen.idx[seq(1, length(cen.idx), by=2)]
  qcen.idx <- cen.idx[seq(2, length(cen.idx), by=2)]
  elementMetadata(cb)$CEN <- 'arm'
  elementMetadata(cb)[cen.idx,]$CEN <- 'cen'
  
  spread <- c(1:spread)
  for(each.spread in spread){
    p.cen.idx <- setdiff((pcen.idx - each.spread), pcen.idx)
    q.cen.idx <- setdiff((qcen.idx + each.spread), qcen.idx)
    
    # Remove cytobands not on the same chromosomes as their corresponding centromere
    rm.idx <- as.character(seqnames(cb[pcen.idx,])) == as.character(seqnames(cb[p.cen.idx,]))
    if(any(!rm.idx)) p.cen.idx <- p.cen.idx[-which(!rm.idx)]
    rm.idx <- as.character(seqnames(cb[qcen.idx,])) == as.character(seqnames(cb[q.cen.idx,]))
    if(any(!rm.idx)) q.cen.idx <- q.cen.idx[-which(!rm.idx)]
    
    elementMetadata(cb[p.cen.idx,])$CEN <- paste0("pcen", each.spread)
    elementMetadata(cb[q.cen.idx,])$CEN <- paste0("qcen", each.spread)
  }
  
  roi <- c(paste0("pcen", rev(spread)), "cen", paste0("qcen", spread))
  list("cb"=cb, "roi"=roi)
}
aggregateGr <- function(list.gr){
  list.gr <- as.list(list.gr)
  # Loop to combine the first two GRanges object of the list, save it to the first element of the list, pop out the second element
  while(length(list.gr) >= 2){
    x <- list.gr[[1]]
    y <- list.gr[[2]]
    # Make a composite GRanges object that contains all possible segments
    int.gr <- disjoin(sort(c(granges(x), granges(y))))
    int.gr <- sort(c(int.gr, gaps(int.gr)))
    
    # Function to create an nrow(int.gr) matrix containing all the mapped elements from x and y
    xy.list <- lapply(list(x, y), function(z, int.gr){
      mat.blank <- matrix(nrow=length(int.gr),
                          ncol=1)
      olaps <- findOverlaps(int.gr, z, type = 'within')
      
      mat.fill <- apply(elementMetadata(z), 2, function(meta.z){
        mat.blank[queryHits(olaps), ] <- meta.z[subjectHits(olaps)]
        mat.blank
      })
      as.data.frame(mat.fill, check.names=FALSE)
    }, int.gr=int.gr)
    
    # Combining the x-y elements back into the intersected GRanges object
    elementMetadata(int.gr) <- cbind(elementMetadata(int.gr), do.call("cbind", xy.list))
    
    list.gr[[2]] <- NULL
    list.gr[[1]] <- int.gr
  }
  list.gr[[1]]
}

# Generates heatmap of peak diff
drawCens <- function(bp.idx, max.y, 
                     cencol='cadetblue3', pcencol='cadetblue1'){
  mid.row <- (nrow(bp.idx)+1)/2
  rect(xleft = bp.idx[mid.row+1, ]-0.5, ybottom = -100,
       xright = bp.idx[mid.row+1, ]+0.5, ytop = 100, 
       col=alpha(pcencol, 0.2), border=NA)
  rect(xleft = bp.idx[mid.row, ]-0.5, ybottom = -100,
       xright = bp.idx[mid.row, ]+0.5, ytop = 100, 
       col=alpha(cencol, 0.3), border=NA)
  rect(xleft = bp.idx[mid.row-1, ]-0.5, ybottom = -100,
       xright = bp.idx[mid.row-1, ]+0.5, ytop = 100, 
       col=alpha(pcencol, 0.2), border=NA)
}
rpkm <- function(seg.gr, treads){
  peak.width <- sum(width(seg.gr)) / 1000
  reads <- sum(rep(seg.gr$peak, width(seg.gr)))
  
  scaling.factor <- (treads / 1000000)
  round((reads / scaling.factor) / peak.width,1)
}  
lohChr <- function(){
  loh=paste0("chr", c(1, 2, 3, 6, 8, 10, 11, 16, 18, 21, 22))
  ret=paste0("chr", c(4, 5, 7, 9, 12, 13, 14, 15, 17, 19, 20))
  
  chrcols <- rep(NA, length(chrs))
  names(chrcols) <- chrs
  chrcols[match(loh, chrs)] <- 'black'
  chrcols[match(ret, chrs)] <- 'white'
  chrcols
}
missegregateChr <- function(chrs=paste0("chr", 1:22)){
  ## Transcribed and estimated using Inkscape Scaling
  ## directly from Figure 1.I of the Worrell Paper
  ## representing % of Cells for single-cell nocodazole treatment
  chrcols <- c(4.75, 2.8, 2.8, 0.8, 0.8, 2.1,
               0, 2.1, 0.8, 1.5, 2.1, 1.5,
               2.1, 0, 0, 0, 0.8, 2.8, 0.8,
               2.1, 0.8, 0.8, 0.8)
  names(chrcols) <- c(chrs, "chrX")
  
  chrcols <- round(chrcols/max(chrcols), 3)
  chrcols
}
makeBox <- function(in.mat){
  in.mat[TRUE] <- 1
  in.mat
}
plotBox <- function(cl.data,dat, ...){
  require(beeswarm)
  boxplot(cl.data, outline=FALSE, las=1, ...)
  beeswarm.points <- beeswarm(cl.data, do.plot=FALSE)
  points(y~x, beeswarm.points, pch=16, col=alpha("black", 0.6))
  
  p  <- t.test(cl.data[[1]], cl.data[[length(cl.data)]])$p.value
  if(p < 0.05){
    y.val <- switch(dat,
                    "t"=1.8,
                    "D"=0.9)
    segments(x0 = 1, y0 = y.val, 
             x1 = length(cl.data), y1 = y.val)
    text(x= 1+(length(cl.data)-1)/2, 
         y=y.val, labels="***", pos = 3)
  }
  axis(side = 1, at=c(1:length(cl.data)), labels=rep("", length(cl.data)))
  if(dat=="D"){
    axis(side = 2, at=seq(0, 1, by=0.2), labels=rep("", 6))
  } else if(dat=="t"){
    axis(side = 2, at=c(-2:2), labels=rep("", 5))
  }
  
  return(p)
}

## KS test to get enriched ROI and reduce/summarize across replicates
ROIenrichment <- function(gr, roi.chr, stat){
  all.esize <- lapply(colnames(elementMetadata(gr)), function(id){
    ## Q (NULL) = RPKM across all regions and all chromosomes
    Q <- as.numeric(elementMetadata(gr)[,id])
    
    sapply(roi.chr, function(rois){
      x <- sapply(rois, function(i) {
        tryCatch({
          ## P = RPKM on chrX (rois) for region I (i)
          P=elementMetadata(i)[,id]
          
          switch(stat,
                 "p"=ks.test(P, Q, alternative='less')$p.value,
                 "D"=ks.test(P, Q, alternative='less')$statistic)
        }, error=function(e) {NA})
      })
      round(x, 2)
    })
  })
  names(all.esize) <- colnames(elementMetadata(gr))
  all.esize
}
reduceROIenrichment <- function(grps, all.esize, t='fishers'){
  reduced.esize <- lapply(grps, function(grp, t){
    ## Implements the metap "sumlog" function in a Reduce fashion
    fishers.sumlog <- function(esize){
      esize <- lapply(esize, function(i, min.p){
        i[i == 0] <- min.p
        i
      }, min.p = 0.00001)
      log.esize <- lapply(esize, log)
      chisq.esize <- -2 * Reduce("+", log.esize)
      df <- 2 * length(idx)
      round(pchisq(chisq.esize, df, lower.tail = FALSE),4)
    }
    
    idx <- grep(grp, names(all.esize))
    red.mat <- switch(t,
                      "avg"=Reduce('+', all.esize[idx]) / length(idx),
                      "fishers"=fishers.sumlog(all.esize[idx]))
    red.mat[len.mat < 2] <- NA
    red.mat
  }, t=t)  # (t=avg, fishers)
  names(reduced.esize) <- grps
  reduced.esize
}
runROIpipeline <- function(gr, roi.chr, grps){
  all.esize <- ROIenrichment(gr, roi.chr, 'p') # stat: p, D
  reduced.esize <- reduceROIenrichment(grps, all.esize, 'fishers') # t: avg, fishers (combines pval)
  
  dval.esize <- ROIenrichment(gr, roi.chr, 'D') # stat: p, D
  reduced.Dsize <- reduceROIenrichment(grps, dval.esize, 'avg') # t: avg (combines D), fishers
  
  list("p"=reduced.esize,
       "D"=reduced.Dsize)
}


## t-test reduction of either p-values, or t-statistics
reduceTtest <- function(roi.chr, stat='p'){
  t.mat <- sapply(roi.chr, function(rois, stat){
    sapply(rois, function(chr.roi){
      if(stat == 'p'){
        tryCatch({
          round(sumz(chr.roi$pval)$p,3)
        }, error=function(e){NA})
      } else if(stat=='stat'){
        tryCatch({
          round(mean(chr.roi$tstat, na.rm=TRUE),3)
        }, error=function(e){NA})
      } else if(stat=='stat.sd'){
        tryCatch({
          round(sd(chr.roi$tstat, na.rm=TRUE),3)
        }, error=function(e){NA})
      }
    })
  }, stat=stat)
  
  t.mat
}

## EM separate enriched peaks
clusterStats <- function(Dmat, tmat, ord=NULL){
  # Run EM if no order is given
  if(is.null(ord)){
    Dvals <- as.numeric(Dmat)
    tstat.vals <- as.numeric(tmat)
    na.idx <- which(is.na(Dvals))
    
    fit <- Mclust(Dvals[-na.idx])
    ord <- fit$classification
    
    D.split <- split(Dvals[-na.idx], f=ord)
    Tstat.split <- split(tstat.vals[-na.idx], f=ord)
  } else {
    Dvals <- lapply(ord, function(chrs.grp){
      idx <- match(chrs.grp, rownames(Dmat))
      as.numeric(Dmat[idx,])
    })
    tstat.vals <- lapply(ord, function(chrs.grp){
      idx <- match(chrs.grp, rownames(tmat))
      as.numeric(tmat[idx,])
    })
    
    na.idx <- lapply(Dvals, function(i) which(is.na(i)))
    D.split <- lapply(seq_along(na.idx), function(idx){
      Dvals[[idx]][-na.idx[[idx]]]
    })
    Tstat.split <- lapply(seq_along(na.idx), function(idx){
      tstat.vals[[idx]][-na.idx[[idx]]]
    })
  }
  
  
  list("D"=D.split,
       "t"=Tstat.split)
}


#########################
#### Setup ChIP Data ####
col.param1 <- c(-2, 0, 0, 1) * 10
col.param2 <- c(-2, 0, 1, 8) * 10
scaling.factor <- 10000000

p <- '1'
project <- switch(p,
                  "1"="Sturgill_CENPA-DAXX",
                  "2"="Hake_chipseq_daxx")
grps <- c("Control", "DAXX")

PDIR <- paste0("/mnt/work1/users/pughlab/projects/NET-SEQ/external_data/", project, "/wig")
setwd(PDIR)
chrs <- paste0("chr", c(1:22, "X", "Y"))

##############
#### MAIN ####
if(project == 'Sturgill_CENPA-DAXX'){
  readcnt <- read.table("sturgil_reads.csv", sep=",", header=TRUE,
                        stringsAsFactors = FALSE)
  readcnt$grp <- gsub("^si", "", readcnt$Condition)
  readcnt$rep <- gsub("^.*_", "", readcnt$Experiment)
  read.spl <- split(readcnt, f=readcnt$gr)
}


cbl <- getCENidx(hg19.cytobands, spread=1)
cb <- cbl[['cb']]; roi <- cbl[['roi']]

files <- list.files(PDIR, pattern="wig$")
if(file.exists("wig_gr.Rdata")){
  print(paste0("Loading pre-existing ", project, " data structure..."))
  load("wig_gr.Rdata")
} else {
  print(paste0("Generating ", project, " GRanges object..."))
  ## Read in each bedgraph file for each file
  wig.gr <- lapply(files, readInWigs)
  names(wig.gr) <- files
  
  mpfiles <- list.files("../", pattern="merged_replicated")
  mp.gw <- lapply(grps, function(grp){
    read.grp <- read.spl[[grp]]
    
    ## Rename the wig files based on replicate name
    grp.gr <- wig.gr[grep(grp, names(wig.gr), ignore.case=TRUE)]
    names(grp.gr) <- names(grp.gr) %>% 
      gsub("^.*_Rep", "Rep", ., ignore.case=TRUE) %>%
      gsub("_.*", "", .)
    
    ## Load in the predefined merged peaks bed file
    f <- mpfiles[grep(grp, mpfiles, ignore.case=TRUE)]
    merged.peaks <- read.table(file.path("..", f), header=FALSE, sep="\t",
                               stringsAsFactors = FALSE)
    colnames(merged.peaks) <- c("chr", "start", "end", "range")
    gr <- makeGRangesFromDataFrame(merged.peaks)
    
    ## Overlap each of the replicates with the merged peaks
    adj.grp.gr <- lapply(names(grp.gr), function(repl){
      g.gr <- grp.gr[[repl]]
      treads <- read.grp[with(read.grp, rep == repl), 'Total.reads']
      
      ov <- findOverlaps(g.gr, gr)
      split.ov <- split(ov, f=subjectHits(ov))
      max.pk <- sapply(split.ov, function(s.ov) max(g.gr[queryHits(s.ov),]$peak))
      gr$peak <- round((max.pk / treads) * 10000000, 1)
      
      gr$rpkm <- sapply(split.ov, function(s.ov) rpkm(seg.gr=g.gr[queryHits(s.ov),],
                                                      treads=treads))

      gr
    })
    names(adj.grp.gr) <- names(grp.gr) 
    
    ## Aggregate the GR: Max peaks/RPKM
    peaks.gr <- aggregateGr(lapply(adj.grp.gr, function(i) i[,1]))
    peaks.gr <- peaks.gr[-which(is.na(peaks.gr[,1]$peak)),]
    colnames(elementMetadata(peaks.gr)) <- names(adj.grp.gr)
    
    rpkm.gr <- aggregateGr(lapply(adj.grp.gr, function(i) i[,2]))
    rpkm.gr <- rpkm.gr[-which(is.na(rpkm.gr[,1]$rpkm)),]
    colnames(elementMetadata(rpkm.gr)) <- names(adj.grp.gr)
    
    list("peaks"=peaks.gr, "rpkm"=rpkm.gr)
  })
  names(mp.gw) <- grps
  
  .aggregateAndCleanGR <- function(grp){
    grp.gr <- aggregateGr(lapply(mp.gw, function(i) i[[grp]]))
    grp.gr <- grp.gr[-which(apply(elementMetadata(grp.gr), 1, 
                                      function(i) any(is.na(i)))),]
    colnames(elementMetadata(grp.gr)) <- paste0(sort(rep(grps, 4)), "_",
                                                  colnames(elementMetadata(grp.gr)))
    grp.gr
  }
  peaks.gr <- .aggregateAndCleanGR('peaks')
  rpkm.gr <- .aggregateAndCleanGR('rpkm')
  
  save(rpkm.gr, peaks.gr, mp.gw, file="wig_gr.Rdata")
}

#### Map RPKM to ROI ####
control <- switch(p,
                  "1"="Control",
                  "2"="H33")
if(p == '1') gr <- rpkm.gr
control2 <- 'H3Y'
roi.chr <- lapply(roi, overlapGrWithROI, gr=rpkm.gr) ## RPKM combined
ctrl.roi.chr <- lapply(roi, overlapGrWithROI, gr=mp.gw[['Control']][['rpkm']]) ## RPKM control
daxx.roi.chr <- lapply(roi, overlapGrWithROI, gr=mp.gw[['DAXX']][['rpkm']]) ## RPKM daxx
names(daxx.roi.chr) <- names(ctrl.roi.chr) <- names(roi.chr) <- roi

#### ROI Enrichment ####
### Checks to see if the centromeric regions are enriched in peaks (RPKM or Max Height)
# Runs a KS-test for each sample, compares each ROI (e.g. centromere) for each chromosome
# against the total background distribution (i.e. all RPKM across entire genome)
# Then combines D-values as average, or p-values using the Fishers method
len.mat <- sapply(roi.chr, function(i) sapply(i, length))
rpkm.roi <- runROIpipeline(rpkm.gr, roi.chr, grps) # reduced.esize, reduced.Dsize
ctrl.rpkm.roi <- runROIpipeline(mp.gw[['Control']][['rpkm']], ctrl.roi.chr, "*") # reduced.esize, reduced.Dsize
daxx.rpkm.roi <- runROIpipeline(mp.gw[['DAXX']][['rpkm']], daxx.roi.chr, "*") # reduced.esize, reduced.Dsize


#### Collapse t-statistics ####
# Collapses the t-test values between Control and Daxx overlapping peaks
# using Fishers method for p, or average/sd for t-stat
t.mat <- reduceTtest(roi.chr, 'p')
tstat.mat <- reduceTtest(roi.chr, 'stat')
tstat.sd.mat <- reduceTtest(roi.chr, 'stat.sd')

#### EM clusters of enrichment ####
control.mat <- rpkm.roi[['D']][['Control']]
em.pericen.dt <- clusterStats(control.mat, tstat.mat)
em.cen.dt <- clusterStats(control.mat[,'cen',drop=FALSE], 
                      tstat.mat[,'cen',drop=FALSE])

split.chrs <- lapply(split(lohChr(), factor(lohChr())), names)
loh.pericen.dt <- clusterStats(control.mat, tstat.mat,
                      split.chrs)
loh.cen.dt <- clusterStats(control.mat[,'cen',drop=FALSE], 
                      tstat.mat[,'cen',drop=FALSE], split.chrs)


p <- list()
pdf(paste0(project, ".enrichBox.pdf"), width = 9, height = 7)
clus.scrn <- split.screen(c(2,6))
# >> [1,1]: EM clusters (periCEN) - D-bal
screen(clus.scrn[2]); par(mar=c(1, 2, 4.1, 0))
p[['em.pericenD']] <- plotBox(em.pericen.dt[['D']], dat='D', 
                              ylim=c(0,1), ylab="D", xaxt='n', main="periCEN")

# [2,1]: EM clusters (periCEN) - t-stat
screen(clus.scrn[8]); par(mar=c(5.1, 2, 1, 0))
p[['em.pericenT']] <- plotBox(em.pericen.dt[['t']], dat='t', 
                              ylim=c(-2, 2), ylab="t-stat")


# >> [1,2]: EM clusters (CEN) - D-bal
screen(clus.scrn[3]); par(mar=c(1, 2, 4.1, 0))
p[['em.cenD']] <- plotBox(em.cen.dt[['D']], dat='D', ylim=c(0,1), 
                          ylab="", xaxt='n', yaxt='n', main="CEN")

# [2,2]: EM clusters (CEN) - t-stat
screen(clus.scrn[9]); par(mar=c(5.1, 2, 1, 0))
p[['em.cenT']] <- plotBox(em.cen.dt[['t']], dat='t', ylim=c(-2, 2), 
                          yaxt='n', ylab="")


# >> [1,3]: LOH clusters (periCEN) - D-bal
screen(clus.scrn[4]); par(mar=c(1, 2, 4.1, 0))
p[['loh.pericenD']] <- plotBox(loh.pericen.dt[['D']][c(2,1)], dat='D', ylim=c(0,1), 
                               ylab="", xaxt='n', yaxt='n', main="periCEN")

# [2,3]: LOH clusters (periCEN) - t-stat
screen(clus.scrn[10]); par(mar=c(5.1, 2, 1, 0))
p[['loh.pericenT']] <- plotBox(loh.pericen.dt[['t']][c(2,1)], dat='t', ylim=c(-2, 2), 
                               ylab="", yaxt='n',names=c("Het", "LOH"))

# >> [1,4]: LOH clusters (CEN) - D-bal
screen(clus.scrn[5]); par(mar=c(1, 2, 4.1, 0))
p[['loh.cenD']] <- plotBox(loh.cen.dt[['D']][c(2,1)], dat='D', ylim=c(0,1), 
                           ylab="", xaxt='n', yaxt='n', main="CEN")

# [2,4]: LOH clusters (CEN) - t-stat
screen(clus.scrn[11]); par(mar=c(5.1, 2, 1, 0))
p[['loh.cenT']] <- plotBox(loh.cen.dt[['t']][c(2,1)], dat='t', ylim=c(-2, 2), 
                           ylab="", yaxt='n', names=c("Het", "LOH"))
close.screen(all.screens=TRUE)
dev.off()






#### Barplot + Line overlay ####
pdf(paste0(project, ".bars.pdf"), width = 10, height = 3.5)
scrn.idx <- split.screen(matrix(c(0, 1, 0, 0.15,
                           0, 1, 0.15, 0.6,
                           0, 1, 0.6, 1), byrow=TRUE, ncol=4))
xaxis.scrn <- scrn.idx[1]
tstat.scrn <- scrn.idx[3]  # ROI per chrom RPKM significance
barplot.scrn <- scrn.idx[2]

## Setup for plots:
control.bp <- t(reduced.esize[['Control']])
rois.idx <- match(c("pcen1", "cen", "qcen1"), rownames(control.bp))
control.bp <- control.bp[rois.idx,]

bp.idx <- barplot(makeBox(control.bp), beside=TRUE, plot = FALSE)

## Plots the t-test statistic for siRNA DAXX compared to control
screen(tstat.scrn); par(mar=c(0.05, 4.1, 1, 2.1))

t.pvals <- t(t.mat)[rois.idx,]
t.stats <- t(tstat.mat)[rois.idx,]
t.sd <- t(tstat.sd.mat)[rois.idx,]

plot(0, type='n', xlim=c(min(bp.idx) - 0.5, 
                         max(bp.idx) + 0.5), 
     ylim=c(-5,3), xaxt='n', xlab='',
     las=1, ylab='t-stat')
drawCens(bp.idx, 3)
abline(h = 0, lty=5, col="black")

sapply(chrs, function(chr){
  print(chr)
  chr.idx <- match(chr, chrs)
  x.range <- bp.idx[, chr.idx]
  
  segments(x0 = x.range, y0 = t.stats[, chr.idx] - t.sd[, chr.idx], 
           x1 = x.range, y1 = t.stats[, chr.idx] + t.sd[, chr.idx])
  points(x=x.range, y=t.stats[, chr.idx],
         pch=16, col=alpha("black", 0.8))
  if(any(t.pvals[, chr.idx] < 0.05, na.rm=TRUE)){
    sig.idx <- which(t.pvals[, chr.idx] < 0.05)
    y.pos <- t.stats[sig.idx, chr.idx] - t.sd[sig.idx, chr.idx] - 0.5
    text(x = x.range[sig.idx], y=y.pos, labels="*", col="red")

    # points(x=x.range[sig.idx], y=t.stats[sig.idx, chr.idx],
    #        pch=16, col=alpha("red", 1))
  }
})
abline(v = (bp.idx[nrow(bp.idx),] + 1), lty=2, col="grey")


## Plots the control RPKM across chromosomes and ROI
screen(barplot.scrn); par(mar=c(0.05, 4.1, 0.05, 2.1))

peaks.roi <- roi.chr[rois.idx]
control.bp <- t(reduced.esize[['Control']])
rois.idx <- match(c("pcen1", "cen", "qcen1"), rownames(control.bp))
control.bp <- control.bp[rois.idx,]

maxy <- 300000
plot(0, type='n', xlim=c(min(bp.idx) - 0.5, 
                         max(bp.idx) + 0.5), 
     ylim=log2(c(1, maxy)), xaxt='n', xlab='',
     las=1, ylab='Log2(RPKM)')
drawCens(bp.idx, log2(maxy))

sapply(chrs, function(chr){
  sapply(roi[rois.idx], function(r){
    roi.idx <- match(r, rownames(control.bp))
    chr.idx <- match(chr, chrs)
    x.range <- bp.idx[roi.idx, chr.idx]
    
    peaks.roi <- elementMetadata(roi.chr[[r]][[chr]])
    # x.space <- round(1/nrow(peaks.roi), 3)
    
    if(nrow(peaks.roi) > 0) {
      grp <- 'Control'; grp.col <- 'grey'
      grp.idx <- grep(grp, colnames(peaks.roi))

      getQuantile <- function(i, q, log.stat=FALSE){
        qi <- quantile(as.numeric(unlist(i)), q)
        if(log.stat) qi <- log2(qi)
        qi
      }
      hiq <- getQuantile(peaks.roi[,grp.idx], 0.75, TRUE)
      loq <- getQuantile(peaks.roi[,grp.idx], 0.25, TRUE)
      med <- getQuantile(peaks.roi[,grp.idx], 0.5, TRUE)
      
      rect(xleft = x.range - 0.5, ybottom = loq,
           xright = x.range + 0.5, ytop = hiq, 
           border=TRUE, col=grp.col, lwd=0.5)
      segments(x0 = x.range - 0.5, y0 = med,
               x1 = x.range + 0.5, y1 = med)
      
      if(is.na(control.bp[roi.idx, chr.idx])) control.bp[roi.idx, chr.idx] <- 1
      if(control.bp[roi.idx, chr.idx] < 0.05){
        text(x = x.range, y=hiq + 1, labels = "*", col="red")
      }
    }
  })
})
abline(v = (bp.idx[nrow(bp.idx),] + 1), lty=2, col="grey")




## Labels chromosomes
screen(xaxis.scrn); par(mar=c(1.1, 4.1, 0.05, 2.1))
plot(0, type='n', xaxt='n', yaxt='n', ylab='', xlab='', axes=FALSE,
     xlim=c(min(bp.idx) - 0.5, max(bp.idx) + 0.5), ylim=c(0,1))
text(x=bp.idx[2,], y=rep(0.5, ncol(bp.idx)), 
     labels = chrs, srt=90, cex=0.7)

close.screen(all.screens=TRUE)
dev.off()




