getChromSnpNeg <- function(x, bin.metric='af.quant.mu'){
  density.af.quant <- density(rep(x[,bin.metric], x[,'num.of.bins']), na.rm=TRUE)
  density.peak <- which(density.af.quant$y == max(density.af.quant$y))[1]
  return(density.af.quant$x[density.peak])
#   weighted.af <- apply(x, 1, function(chrom.row){
#     (as.integer(as.character(chrom.row['num.of.bins'])) * 
#        as.numeric(as.character(chrom.row['af.quant.mu'])))
#   })
#   chr.weighted.af <- sum(unlist(weighted.af), na.rm=TRUE) / sum(x[,'num.of.bins'], na.rm=TRUE)
#   return(chr.weighted.af)
}

getChromCnNeg <- function(x){
  density.af.quant <- density(rep(x[,'copy.number'], x[,'num.of.bins']), na.rm=TRUE)
  density.peak <- which(density.af.quant$y == max(density.af.quant$y))[1]
  return(density.af.quant$x[density.peak])
}

getMeanAndSd <- function(x){
  mean.val <- mean(x, na.rm=TRUE)
  sd.val <- sd(x, na.rm=TRUE)
  return(paste(round(mean.val, 2), "+/-", round(sd.val, 2), sep=" "))
}

fitCnLabels <- function(fit.df, diploid.stat=TRUE){
  # Find the highest peak and the zero-peak for reference
  average.sd <- mean(fit.df[which(fit.df$mu > 100),]$sigma)
  diploid.peak <- fit.df[which(fit.df$pi == max(fit.df$pi)),]$mu
  #null.peak <- fit.df[1,]$mu
  null.peak <- 0
  if(null.peak > 50){
    print(paste("Warning: Null peak is greater than 50 - ", null.peak, sep=""))
  }
  
  # If diploid assumption, assume highest peak and null difference is 2*
  if(diploid.stat){
    peak.diff <- (diploid.peak - null.peak)/2
  } else {
    peak.diff <- (diploid.peak - null.peak)/1
  }
  
  # Estimate CN and measure difference between estimate and observed
  rounded.cn <- round(fit.df$mu / peak.diff, 0)
  rounded.cn[which(rounded.cn < 0 )] <- 0
  estimate.diff <- apply(data.frame(pred.val=(rounded.cn * peak.diff),
                                    raw.val=fit.df$mu), 1, diff)
  # Obtain the -log probability of obtaining just that CN value
  prob.model <- -dnorm(x = (rounded.cn * peak.diff),
                       mean=fit.df$mu,
                       sd=fit.df$sigma, log=TRUE)
  
  estimated.cn <- cbind(fit.df,
                        abs.cn=rounded.cn,
                        diff=estimate.diff,
                        log.prob=prob.model)
  return(estimated.cn)
}

#returns the data.frame needed for aggregateLOH
genDiseaseDF <- function(x, type='bin'){
  if(type %in% 'bin'){
    quant.type <- 'quantile'
  } else if (type %in% 'seg'){
    quant.type <- bin.m
  } else if (type %in% 'nosnp'){
    quant.type <- NA
  }
  
  chromDisease.list <- t(apply(x, 1, function(x.row){
    if(is.na(x.row[quant.type])){
      loh.stat <- 0
    } else if(x.row[quant.type] <= ecdf.thresh) {
      loh.stat <- 1 
    } else {
      loh.stat <- 0
    }
    
    ret.val <- c(Chromosome = as.character(x.row['chrom']),
                 Start.bp = as.integer(x.row['start.pos']),
                 End.bp = as.integer(x.row['end.pos']),
                 modal_A1 = as.integer(x.row['copy.number']),
                 modal_A2 = 0,
                 LOH=loh.stat)
    return(ret.val)
  }))
  
  return(chromDisease.list)
}

#collapses all segments that have the same copy-state and LOH status to avoid plotting a ton of segments
collapseCopyState <- function(x){
  x <- as.data.frame(x)
  pasted.order <-  paste(x$Chromosome, x$modal_A1, x$LOH, sep="_")
  end.vals <- which(diff(as.numeric(interaction(pasted.order))) != 0)
  start.vals <- c(1, (end.vals + 1))
  end.vals <- c(end.vals, dim(x)[1])
  
  collapse.x <- cbind(x[start.vals, c('Chromosome', 'Start.bp')],
                      x[end.vals, 'End.bp'],
                      x[start.vals, c('modal_A1', 'modal_A2', 'LOH')])
  colnames(collapse.x) <- colnames(x)
  return(as.data.frame(collapse.x))
}

# Uses rle() to obtain consecutive values in a dataframe and index them
getRleIdx <- function(x, col.id=NA, na.val=-100){
  if(is.vector(x)){
    reformat.na.x <- as.character(x)
    
  } else if(is.data.frame(x)){
    # Handles multiple columns
    if(length(col.id) > 1) {
      uniq.id <- apply(x, 1, function(y) paste(y[col.id], collapse="-"))
      x$uniq <- uniq.id
      col.id <- 'uniq'
    }
    reformat.na.x <- as.character(x[,col.id])
  }
  
  # Fill in NA values
  reformat.na.x[which(is.na(reformat.na.x))] <- na.val
  rle.x <- rle(reformat.na.x)
  
  #Get the array index for the start-to-end of each unique value/changepoint
  rle.x$start.idx <- c(1, (cumsum(rle.x$lengths) + 1)[-length(rle.x$lengths)])
  rle.x$end.idx <- rle.x$start.idx + (rle.x$lengths - 1)
  rle.x$values[which(rle.x$values == na.val)] <- NA
  rle.x$na.stat <- !is.na(rle.x$values) 
  
  return(rle.x)
}

# Uses rle() to obtain consecutive values in a dataframe and index them
getBlankVcfCol <- function(x, col.id=NA, na.val=-100){
  
}