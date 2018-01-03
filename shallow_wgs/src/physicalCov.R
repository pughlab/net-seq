getEMestimates <- function(wgsCov.df, break.num=1600, est.mu = NA, est.sig = NA){
  # Obtain the histogram and number of counts within each histogram bin
  segDoc.his <-  hist(wgsCov.df$count, breaks=break.num, xlim = c(0,5000))
  segDoc.df <- data.frame(mid=segDoc.his$mids, cou=segDoc.his$counts)  
  
  # Estimate the number of mu's and sigma's based on KDE
  if(all(c(is.na(est.mu), is.na(est.sig)))){
    # KDE of genome counts to get rough estimate of local maximas
    segDoc.dens <- density(wgsCov.df$count, na.rm=TRUE, bw='SJ', kernel='gaussian', n=512) 
    seg.min.max <- getLocalMinMax(segDoc.dens, 'lower', 0.01)
    
    # Estimate parameters for the mixed gaussian distribution
    estMean <- segDoc.dens$x[seg.min.max$max]
    estSigma <- getSigmaEstimate(estMean, segDoc.dens$x[seg.min.max$min], max.val=300000)
  } else {
    # Take input as rough approximation
    estMean <- est.mu
    estSigma <- est.sig
  }
  
  segDoc.df <- segDoc.df[which(segDoc.df$mid < (estMean[length(estMean)] + (4 * estSigma[length(estSigma)]))),]
  estDist <- 'norm'
  fitpro <- mix(as.mixdata(segDoc.df), mixparam(mu=estMean, sigma=estSigma), dist=estDist)
#   estDist <- 'binom'
#   fitpro <- mix(as.mixdata(segDoc.df), mixparam(mu=estMean, sigma=estSigma), 
#                 dist=estDist, constr=mixconstr(consigma='BINOM', size=rep(500000, length(estMean))))
  return(fitpro)
}



getCnSegments <- function(wgsCov.df){
  # Collapses all 500k bins into segments based on no change between consecutive bins
  switch.points <- which(diff(as.integer(wgsCov.df$copy.number)) != 0)
  seg.bin.length <- diff(switch.points)
  
  collapse.seg.list <- lapply(c(0:length(switch.points)), function(x) {
    collapseSegments(x, switch.points, seg.bin.length, wgsCov.df)
  })

  
  collapse.seg.df <- do.call("rbind", collapse.seg.list)
}

collapseSegments <- function(x, switch.points, seg.bin.length, wgsCov.df, bin.size=500000){
  if(x == 0 && (length(switch.points) > 0)){
    # Start from row 1 until first switch.point
    srow <- 1
    erow <- switch.points[1]
  } else if (length(switch.points) == 0){
    srow <- 1
    erow <- nrow(wgsCov.df)
  } else if (x == length(switch.points)){
    # Start from last switch.point until last row of dataframe
    srow <- (switch.points[x] + 1)
    erow <- dim(wgsCov.df)[1]
  } else {
    # Start from switch.point to the next switch.point based on distance between
    srow <- (switch.points[x] + 1)
    erow <- (switch.points[x] + seg.bin.length[x])
  }
  #Assemble seg statistics
  startSeg <- wgsCov.df[srow, c('chrom', 'start.pos')]
  endSeg <- wgsCov.df[erow, c('end.pos', 'copy.number')]
  count.df <- data.frame(num.of.bins=round((endSeg$end.pos - startSeg$start.pos) / bin.size, 0),
                         count=mean(wgsCov.df[c(srow:erow), 'count']),
                         sd.count=sd(wgsCov.df[c(srow:erow), 'count']),
                         avg.prob=mean(wgsCov.df[c(srow:erow), as.character(endSeg$copy.number)]),
                         sd.prob=sd(wgsCov.df[c(srow:erow), as.character(endSeg$copy.number)]))
  if(!is.na(match("quantile", colnames(wgsCov.df)))){
    count.df <- cbind(count.df,
                      data.frame(af.quant.mu=mean(wgsCov.df[c(srow:erow), 'quantile'], na.rm = TRUE),
                                 af.quant.med=median(wgsCov.df[c(srow:erow), 'quantile'], na.rm=TRUE)))
  }
  seg.df <- cbind(startSeg, endSeg, count.df)
  return(seg.df)
}

collapseCov <- function(wgsCov.df, seg.limit=NA){
  wgsCov.list <- split(wgsCov.df, wgsCov.df$chrom)
  chr.seg.list <- lapply(wgsCov.list, function(wgsCov.df) getCnSegments(wgsCov.df))
  collapse.seg.df <- do.call("rbind", chr.seg.list)
  
  if(!is.na(seg.limit)){
    collapse.seg.df <- collapse.seg.df[which(collapse.seg.df$num.of.bins > seg.limit),]
    collapse.seg.df <- collapseCov(collapse.seg.df)
  }
  
  return(collapse.seg.df)
}

getMaxProb <- function(x, column.names){
  column.names[which(x == max(x))]
}


dbinomCounts <- function(x, success.cnt, bsize=500000){
  success.prob <- (success.cnt / bsize)
  dbinom(as.integer(as.character(x['count'])), bsize, success.prob)
}

