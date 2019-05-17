library(GenomicRanges)
library(PLTK)
library(dplyr)
library(RColorBrewer)
library(scales)
library(metap)
library(mclust)
require(matrixStats)

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
    names(Dvals) <- names(tstat.vals) <- gsub("\\.D.*$", "", rownames(Dmat))
    
    set.seed(1234)
    fit <- Mclust(Dvals[-na.idx])
    ord <- fit$classification
    
    D.split <- split(Dvals[-na.idx], f=ord)
    Tstat.split <- split(tstat.vals[-na.idx], f=ord)
  } else {
    Dvals <- lapply(ord, function(chrs.grp){
      idx <- which(rownames(Dmat) %in% chrs.grp)
      Dmat[idx,]
    })
    tstat.vals <- lapply(ord, function(chrs.grp){
      idx <- which(rownames(tmat) %in% chrs.grp)
      tmat[idx,]
    })

    na.idx <- lapply(Dvals, function(i) which(is.na(i)))
    D.split <- lapply(seq_along(na.idx), function(idx){
      Dvals[[idx]][-na.idx[[idx]]]
    })
    Tstat.split <- lapply(seq_along(na.idx), function(idx){
      tstat.vals[[idx]][-na.idx[[idx]]]
    })
    names(D.split) <- names(Tstat.split) <- names(ord)
  }
  
  
  list("D"=D.split,
       "t"=Tstat.split)
}
clusterEnrichment <- function(control.mat, tstat.mat, r='cen', 
                              ord=NULL, cols=c("#a6cee3", "#1f78b4", "#b2df8a")){
  em.cen.dt <- clusterStats(control.mat[,r,drop=FALSE], 
                            tstat.mat[,r,drop=FALSE], ord)
  em.cen.dt <- data.frame("D"=unlist(em.cen.dt[['D']]),
                          't'=unlist(em.cen.dt[['t']]),
                          'grp'=rep(names(em.cen.dt[['D']]), 
                                    sapply(em.cen.dt[['D']], length)))
  em.cen.dt$cols <- cols[as.integer(factor(em.cen.dt$grp))]
  em.cen.dt$pch <- 16
  em.cen.dt$chr <- gsub("^.*\\.", "", rownames(em.cen.dt))
  em.cen.dt
}

## Scatterplot with correlations
plotScatPlotCor <- function(comp.df, add=FALSE, text.y=NULL, 
                            text.x=1, id=NULL, ...){
  ct <- cor.test(comp.df[,1], comp.df[,2])
  if(is.null(text.y)) text.y <- max(comp.df[,2], na.rm=TRUE) * 0.99
  
  if(!add){
    plot(comp.df, ...)
    if(!is.null(id)) id <- paste0(id, ": ")
    text(text.x, text.y, pos=2,
         labels = paste0(id, "r=", round(ct$estimate,2),
                         ", p=", round(ct$p.value,2)), ...)
  } else {
    points(comp.df, ...)
    if(!is.null(id)) id <- paste0(id, ": ")
    text(text.x, text.y,
         labels = paste0(id, "r=", round(ct$estimate,2),
                         ", p=", round(ct$p.value,2)),pos=2, ...)
  }
  text.y
}
## RPKM Plots
# Gets the relative size of each centromere + surrounding regions
# and returns a scaled down granges with a set bin.size
getRelativeRoiSize <- function(gr.list, chrs, min.size){
  chr.sizes <- as.data.frame(t(sapply(chrs, function(chr){
    se.pos <- sapply(seq_along(gr.list), function(gr.idx){
      c(min(start(gr.list[[gr.idx]][[chr]])), 
        max(end(gr.list[[gr.idx]][[chr]])))
    })
    se.pos[is.infinite(se.pos)] <- NA
    c("start"=min(se.pos), "end"=max(se.pos), 
      "size"=max(se.pos) - min(se.pos))
  })))
  
  chr.sizes[which(chr.sizes$size < min.size),]$size <- min.size
  rel.size <- with(chr.sizes, ceiling(1/(min(size, na.rm=TRUE) / size)))
  chr.sizes$rsize <- rel.size 
  chr.sizes$chr <- rownames(chr.sizes)
  chr.sizes[is.na(chr.sizes)] <- min(chr.sizes$size, na.rm=TRUE)
  chr.sizes$csize <- c(1, cumsum(chr.sizes$size)[-nrow(chr.sizes)])
  
  chr.bins <- apply(chr.sizes, 1, function(chr.row){
    .ai <-function(x) {as.integer(as.character(x))}
    binsize <- (.ai(chr.row['end']) - .ai(chr.row['start'])) / .ai(chr.row['rsize'])
    pos <- seq(.ai(chr.row['start']), .ai(chr.row['end']), by=binsize)
    if(binsize == 0){
      makeGRangesFromDataFrame(data.frame('chrom'=as.character(chr.row['chr']), 
                                          'start'=1, 'end'=min.size, 'size'=1),
                               keep.extra.columns = TRUE)
    } else {
      makeGRangesFromDataFrame(data.frame('chrom'=rep(as.character(chr.row['chr']), 
                                                      length(pos)-1), 
                                          'start'=pos[-length(pos)], 
                                          'end'=pos[-1]+1,
                                          'size'=rep(as.character(chr.row['rsize']), 
                                                     length(pos)-1)),
                               keep.extra.columns = TRUE)
    }
  })
  list("size"=chr.sizes, "bin"=chr.bins)
}
# Aggregates the GRanges for the reference bins and the ROI GRanges
aggGr <- function(chr, gr1, gr2){
  ov <- findOverlaps(gr1[[chr]], gr2[['bin']][[chr]])
  emgr1 <- elementMetadata(gr1[[chr]])
  emgr2 <- elementMetadata(gr2[['bin']][[chr]])
  
  spl.meta <- split(as.data.frame(emgr1[queryHits(ov),]), 
                    f=subjectHits(ov))
  collapse.meta <- lapply(spl.meta, function(i) colMeans2(as.matrix(i)))
  
  emgr2 <- matrix(nrow=nrow(emgr2), 
                  ncol=ncol(emgr1))
  for(idx in seq_along(names(spl.meta))){
    ridx <- as.integer(names(spl.meta)[idx])
    emgr2[ridx,] <- collapse.meta[[idx]]
  }
  colnames(emgr2) <- colnames(emgr1)
  elementMetadata(gr2[['bin']][[chr]]) <- emgr2
  gr2[['bin']][[chr]]
}
# Plots
plotChrChip <- function(chrs, bin.gr, cols, grp){
  .drawCen <- function(chr, minstart, maxend, spacer){
    cen.chr <- cb[which(seqnames(cb) == chr & cb$CEN == 'cen'),]
    cen.loc <- c(min(start(cen.chr)), max(end(cen.chr)))
    
    cen.st <- cen.loc[1] - minstart + spacer
    cen.end <- cen.loc[2] - minstart + spacer
    min.s <- minstart - minstart + spacer
    max.e <- maxend - minstart + spacer
    rect(xleft = if(cen.st > min.s) cen.st else min.s, ybottom = -100, 
         xright = if(cen.end > max.e) max.e else cen.end, ytop = 100,
         col=alpha("grey", 0.2), border=NA)
  }
  
  plot(0, type='n', xaxt='n', xlab='', ylab='log10(RPKM)',
       xlim=c(min(bin.size$csize), max(bin.size$csize)), 
       ylim=if(length(grp) == 1) c(0,6) else c(-4,4))
  axis(side = 1, at = bin.size$csize + (bin.size$size / 2), 
       labels=bin.size$chr, tick = FALSE, las=2, cex.axis=0.7)
  axis(side = 1, at = bin.size$csize, labels=rep("", nrow(bin.size)))
  
  for(chr in chrs){
    chr.bins <- as.data.frame(bin.gr[[chr]])
    chr.spacer <- bin.size[match(chr, bin.size$chr),]$csize
    
    if(length(grp) == 2){
      ## Run a comparison between groups
      col1.idx <- grep(grp[1], colnames(chr.bins))
      col2.idx <- grep(grp[2], colnames(chr.bins))
      expr.val <- rowMeans(chr.bins[,col1.idx]) - rowMeans(chr.bins[,col2.idx])
      col.idx <- c(col1.idx, col2.idx)
      abline(h = 0, lty=2)
    } else {
      col.idx <- grep(grp[1], colnames(chr.bins))
      expr.val <- rowMeans(chr.bins[,col.idx])
    }
    
    if(any(!is.na(chr.bins[,col.idx]))){
      bin.s <- chr.bins$start
      bin.e <- chr.bins$end
      .drawCen(chr, min(bin.s), max(bin.e), chr.spacer)
      
      expr.sign <- integer(length = length(expr.val))
      expr.sign[which(expr.val > 0)] <- 1
      expr.sign[which(expr.val < 0)] <- -1
      rpkm.df <- data.frame(x1=(bin.s - min(bin.s) + chr.spacer), 
                            x2=(bin.e - min(bin.s) + chr.spacer), 
                            y1=0, 
                            y2=expr.sign * log10(abs(expr.val)))
      with(rpkm.df, rect(xleft = x1, ybottom = y1, 
                         xright = x2, ytop = y2, 
                         border=NA, col=cols[1]))        
    }
  }
}

#########################
#### Setup ChIP Data ####
daxx.col <- '#e7298a'
control.col <- '#a6761d'

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
# split chromosomes based on LOH or Het
split.chrs <- lapply(split(lohChr(), factor(lohChr())), names)
names(split.chrs) <- c("LOH", "Het")
loh.cols <- c("#33a02c", "#fb9a99")
loh.dt <- do.call("rbind", lapply(roi, function(r){
  tmp <- clusterEnrichment(control.mat, tstat.mat, r=r,
                    ord=split.chrs, cols=loh.cols)
  if(r!='cen') tmp$pch <- 1
  tmp
}))

# Visualization
pdf(paste0(project, ".enrichBox.pdf"), width = 5, height = 5)
plot(D~t, data=loh.dt, col=loh.dt$cols, pch=loh.dt$pch, 
     xlim=c(-1.5, 1.5), ylim=c(0,0.7), main="LOH/Het chr", 
     ylab="CENPA enrichment (D)", xlab="DAXX-KO CENPA change (t)")
loh.cen.dt <- loh.dt[which(loh.dt$pch == 16),]
text(D~t, data=loh.cen.dt, labels=gsub("chr", "", loh.cen.dt$chr), pos=4)
abline(v = 0, lty=2, col="grey")
legend("bottomleft", legend = levels(loh.dt$grp), fill=rev(loh.cols), horiz=TRUE)
dev.off()



#### Scatterplots: Correlation of features to Missegregation ####
## Correlation of number of CENPA peaks with missegregation rate (CONTROL then DAXX)
ctrl.cnt <- sapply(ctrl.roi.chr, function(r) sapply(r, length))
daxx.cnt <- sapply(daxx.roi.chr, function(r) sapply(r, length))
cenpa.diff <- rowDiffs(cbind(ctrl.cnt[,idx], daxx.cnt[,idx]))

pdf(paste0(project, ".numPeaks.scatterCor.pdf"), width = 5, height = 5)
idx='cen'
comp.df <- data.frame(missegregateChr(), daxx.cnt[1:23,idx])
text.y <- plotScatPlotCor(comp.df, add=FALSE, col=alpha(control.col, 0.7), 
                          pch=16, id='siRNA Control',
                          xlab="Misseg. Fraction", ylab="Number of CENPA peaks")
comp.df <- data.frame(missegregateChr(), ctrl.cnt[1:23,idx])
plotScatPlotCor(comp.df, add=TRUE, col=alpha(daxx.col, 0.7), 
                pch=16, text.y=text.y * 0.95, id='siRNA DAXX')
comp.df <- data.frame(missegregateChr(), cenpa.diff[1:23])
plotScatPlotCor(comp.df, add=TRUE, col=alpha("black", 0.7), pch=17, 
                text.y=text.y * 0.9, id='DAXX - Control')
dev.off()

## Correlation of centromere size with missegregation rate
pdf(paste0(project, ".cenSize.scatterCor.pdf"), width = 5, height = 5)
cen.len <- sapply(roi, function(r) {
  cb.tmp <- cb[which(cb$CEN == r)]
  sapply(split(cb.tmp, f=seqnames(cb.tmp)), function(i){
    sum(width(i))
  })
})
text.y <- plotScatPlotCor(data.frame(missegregateChr(), cen.len[1:23,'cen']), 
                          add=FALSE, col='black', pch=16, id='CEN',
                          xlab="Misseg. Fraction", ylab="Centromere size")
plotScatPlotCor(data.frame(missegregateChr(), cen.len[1:23,'pcen1']),
                add=TRUE, col='red', pch=1, text.y=text.y * 0.97, id='periCEN (p-arm)')
plotScatPlotCor(data.frame(missegregateChr(), cen.len[1:23,'qcen1']),
                add=TRUE, col='blue', pch=1, text.y=text.y * 0.94, id='periCEN (q-arm)')
dev.off()

## Correlation of CENPA expression to missegregation rate
.averageRpkm <-function(rc){ log2(mean(rowMeans(as.matrix(elementMetadata(rc)))))}
daxx.rpkm <- sapply(daxx.roi.chr, function(r) sapply(r, .averageRpkm))
ctrl.rpkm <- sapply(ctrl.roi.chr, function(r) sapply(r, .averageRpkm))
daxx.ctrl.rpkm <- sapply(roi.chr, function(r) sapply(r, .averageRpkm))

pdf(paste0(project, ".cenpaRPKM.scatterCor.pdf"), width = 5, height = 5)
text.y <- plotScatPlotCor(data.frame(missegregateChr(), daxx.rpkm[1:23,'cen']), 
                          add=FALSE, col=alpha(daxx.col, 0.7), 
                          pch=16, id='siRNA DAXX',
                          xlab="Misseg. Fraction", ylab="Average RPKM")
plotScatPlotCor(data.frame(missegregateChr(), ctrl.rpkm[1:23,'cen']),
                add=TRUE, col=alpha(control.col, 0.7), pch=16, 
                text.y=text.y * 0.965, id='siRNA control')
plotScatPlotCor(data.frame(missegregateChr(), daxx.ctrl.rpkm[1:23,'cen']),
                add=TRUE, col=alpha("black", 0.7), pch=17, 
                text.y=text.y * 0.93, id='overlap peaks')
dev.off()



#### Barplot + Line overlay ####
.combineROI <- function(chr, gr){
  c(gr[[1]][[chr]], gr[[2]][[chr]], gr[[3]][[chr]])
}

# Combine the CEN with the periCEN regions
chrs <- names(roi.chr[[1]])
all.roi.chr <- sapply(chrs, .combineROI, gr=roi.chr)
all.ctrl.roi.chr <- sapply(chrs, .combineROI, gr=ctrl.roi.chr)
all.daxx.roi.chr <- sapply(chrs, .combineROI, gr=daxx.roi.chr)

# Create a reference bin set of a set bin-size
gr.list <- list(all.roi.chr, all.ctrl.roi.chr, all.daxx.roi.chr)
gr.ref.bins <- getRelativeRoiSize(gr.list, chrs, min.size=500000)

# Overlap RPKM gr with the reference bins
bin.roi.chr <- sapply(chrs, aggGr, gr1=all.roi.chr, gr2=gr.ref.bins)
bin.ctrl.chr <- sapply(chrs, aggGr, gr1=all.ctrl.roi.chr, gr2=gr.ref.bins)
bin.daxx.chr <- sapply(chrs, aggGr, gr1=all.daxx.roi.chr, gr2=gr.ref.bins)

bin.size <- gr.ref.bins[['size']]

pdf(paste0(project, ".bars.pdf"), width = 10, height = 3.5)
plotChrChip(chrs=chrs, bin.gr=bin.daxx.chr, 
            col=daxx.col, grp='Rep')
plotChrChip(chrs=chrs, bin.gr=bin.ctrl.chr, 
            col=control.col, grp='Rep')
plotChrChip(chrs=chrs, bin.gr=bin.roi.chr, 
            col="black", 
            grp=c('DAXX', 'Control'))
dev.off()











