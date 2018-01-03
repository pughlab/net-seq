compareTumorNormal <- function(seg.list.x, seg.list.y, norm.process = 'median', plot.dens=FALSE, change.perc=NA){  
  # Get the mean of cbs seg A and cbs seg B on a 0-0.5 AF scale, then normalize
  seg.x.df <- getSegMean(seg.list.x)
  seg.y.df <- getSegMean(seg.list.y)
  
  # Closes the gap between tumor and normal means by a given percentage of each bin distance
  if(!is.na(change.perc)){
    if(dim(seg.y.df)[1] == dim(seg.x.df)[1]){
      change.v <- sapply(c(1:dim(seg.x.df)[1]), function(x, change.perc) {
        change.vector <- -((seg.x.df[x,'seg.ab'] - seg.y.df[x,'seg.ab']) * change.perc)
      }, change.perc=change.perc)
      seg.x.df$seg.ab <- seg.x.df$seg.ab + change.v
    } else {
      print("Warning: seg.x.df and seg.y.df are not the same length")
    }
  }
  
  
  seg.xy.dens <- normalizeDistributions(seg.x.df$seg.ab, seg.y.df$seg.ab, norm.process)
  seg.x.df$seg.ab.norm <- seg.xy.dens[['seg.x']]
  seg.y.df$seg.ab.norm <- seg.xy.dens[['seg.y']]
  if(plot.dens){
    plotTNnormalization(seg.x.df, seg.y.df)
  }
  
  # Get the mean sd for each bin
  seg.x.df$ab.sd <- getSdMean(seg.list.x)
  seg.y.df$ab.sd <- getSdMean(seg.list.y)
  
  #Overlap the two seg.x and seg.y
  seg.x.df$id <- paste(seg.x.df$chr, seg.x.df$start.pos, sep="-")
  seg.y.df$id <- paste(seg.y.df$chr, seg.y.df$start.pos, sep="-")
  seg.x.df <- seg.x.df[which(seg.x.df$id %in% seg.y.df$id),]
  seg.y.df <- seg.y.df[which(seg.y.df$id %in% seg.x.df$id),]
  
  #Tests for statistical differences between each bin
  bin.p.val <- sapply(c(1:dim(seg.x.df)[1]), function(x) {
    raw.t.test(c("x"=seg.x.df[x, 'seg.ab.norm'], "sd"=seg.x.df[x, 'ab.sd'], "n"=seg.x.df[x, 'n']), 
               c("x"=seg.y.df[x, 'seg.ab.norm'], "sd"=seg.y.df[x, 'ab.sd'], "n"=seg.y.df[x, 'n']))
    })
  bin.p.val.adj <- p.adjust(bin.p.val, method="fdr", 
                            n=length(bin.p.val[which(! is.na(bin.p.val))]))
  
  # Assemble the dataframe to return
  bin.pval.df <- cbind(seg.x.df[,c("chr", "arm", "start.pos", "end.pos", "n", "seg.ab.norm", "ab.sd")],
                       seg.y.df[,c("n", "seg.ab.norm", "ab.sd")],
                       bin.p.val, bin.p.val.adj)
  colnames(bin.pval.df) <- c("chr", "arm", "start.pos", "end.pos", "n.T", "mean.seg.T", "mean.sd.T",
                             "n.N", "mean.seg.N", "mean.sd.N", "pval", "pval.adj")
  
  return(bin.pval.df)
}  

getSdMean <- function(seg.list){
  #Generate average SD value for A and B allele
  seg.df <- do.call("rbind", seg.list)

  a.sd <- data.frame(abs(seg.df$low.sd.a - seg.df$mean.a), abs(seg.df$high.sd.a - seg.df$mean.a))
  b.sd <- data.frame(abs(seg.df$low.sd.b - seg.df$mean.b), abs(seg.df$high.sd.b - seg.df$mean.b))
  
  ab.sd.df <- data.frame(apply(a.sd, 1, mean), apply(b.sd, 1, mean))
  ab.sd <- apply(ab.sd.df, 1, mean)
  
  return(ab.sd)
}

getSegMean <- function(seg.list){
  #Generate average seg value for A and B allele; remove -1 fillers with NA for removal
  seg.df <- do.call("rbind", seg.list)
  seg.df[which(seg.df$seg.a < 0), 'seg.a'] <- NA
  seg.df[which(seg.df$seg.b < 0), 'seg.b'] <- NA
  
  abs.sega <- c(0.5 - abs(0.5 - seg.df$seg.a))
  abs.segb <- c(0.5 - abs(0.5 - seg.df$seg.b))
  x.df <- data.frame(abs.sega, abs.segb)
  seg.df$seg.ab <- apply(x.df, 1, mean)
  
  return(seg.df)
}

normalizeDistributions <- function(sx.dist, sy.dist, norm.stat){
  require(plyr)
  cn.mat <- t(rbind.fill.matrix(t(matrix(sx.dist)), t(matrix(sy.dist))))
  
  if(norm.stat %in% 'quantile'){
    print("Quantile normalizing the segments")
    require(preprocessCore)
    
    cn.qmat <- normalize.quantiles(cn.mat, copy=TRUE)
  } else if (norm.stat %in% 'median'){
    require(limma)
    
    print("Median normalizing the segments")
    cn.qmat <- normalizeMedianAbsValues(cn.mat)
  } else if (is.na(norm.stat)){
    print("No normalization being used")
    cn.qmat <- cn.mat
  }
 
  
  return(list("seg.x" = cn.qmat[,1],
              "seg.y" = cn.qmat[,2]))
}

  


