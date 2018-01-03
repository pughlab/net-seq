# Runs the Pipeline that outputs the LOH regions and tally-counts
getLohOverlap <- function(all.net.loh, 
                          group.id1=c('NET-001', 'NET-003', 'NET-008', 'NET-009'), 
                          group.id2=c('NET-002', 'NET-006', 'NET-007')){
  print("Running LOH Overlap pipeline")
  populated.interval <- list()
  for(each.chr in paste("chr", c(1:22, "X"), sep="")){
    all.net.interval <- lapply(all.net.loh, createIntervals, each.chr)
    all.chr.intersect <- intersectAllIntervals(all.net.interval)
    if(nrow(all.chr.intersect) != 0){
      print(paste("Populating: ", each.chr))
      populated.interval[[each.chr]] <- populateIntervals(sample.int=all.net.interval, 
                                                          ref.int=all.chr.intersect, 
                                                          group.stat=TRUE, chr=each.chr, 
                                                          group.id1=group.id1, group.id2=group.id2)
    }
  }
  populated.interval.collapse <- do.call("rbind", populated.interval)
  populated.interval <- populated.interval.collapse
  return(populated.interval)
}

# Creates a matrix of all possible intersections between the intervals
intersectAllIntervals <- function(all.interval){
  all.endpoints <- as.integer(as.character(unique(c(0, unlist(lapply(all.interval, function(x) x[,c(1,2)]))))))
  all.endpoints <- sort(all.endpoints)
  
  all.intervals <- matrix(c(all.endpoints[-length(all.endpoints)],
                            all.endpoints[-1]), ncol=2)
  return(all.intervals)
}

# Creates intervals from each LOH segment across all Samples
createIntervals <- function(x, each.chr){
  x[[each.chr]]$start <- as.integer(as.character(x[[each.chr]]$start))
  x[[each.chr]]$end <- as.integer(as.character(x[[each.chr]]$end))
  cn.state <- x[[each.chr]]$CNt
  chr.interval.net <- tryCatch({
    Intervals(as.matrix(x[[each.chr]][,c("start", "end")]))
  }, error=function(e){Intervals(matrix(c(0,0), ncol=2))})
  chr.interval.net <- cbind(chr.interval.net, x[[each.chr]]$CNt)
  
  return(chr.interval.net)
}

# Creates a dataframe populating the "all-endpoint" matrix with counts of LOH overlap
populateIntervals <- function(sample.int, ref.int, group.stat=FALSE, chr, group.id1, group.id2){
  # Create base dataframe
  tally.df <- as.data.frame(ref.int)
  tally.df <- cbind("chr"=rep(chr, nrow(tally.df)), tally.df)
  tally.df$length <- tally.df$V2 - tally.df$V1
  tally.df[,c(paste('all_(n=', length(group.id1) + length(group.id2), ")", sep=""),
              paste('PNETs_(n=', length(group.id1), ")", sep=""), 
              paste('GINETs_(n=', length(group.id2), ")", sep=""))] <- 0
  colnames(tally.df)[c(2:3)] <- c("start.pos", "end.pos")
  
  # Create a running tally of the LOH segments in each interval
  ref.interval <- Intervals(ref.int)
  for(each.sample in names(sample.int)){
    cnt.stat <- sapply(interval_overlap(ref.interval, 
                                        Intervals(sample.int[[each.sample]][,c(1,2)])), 
                       function(cnt) length(cnt) > 0)
    
    # Tally up the group count
    tally.df[,5] <- tally.df[,5] + cnt.stat
    if(each.sample %in% group.id1 & group.stat){
      tally.df[,6] <- tally.df[,6] + cnt.stat
    } else if(each.sample %in% group.id2 & group.stat){
      tally.df[,7] <- tally.df[,7] + cnt.stat
    } else {
      print(paste("Group tally mode is set to:", group.stat))
      print(paste("Each sample: ", each.sample))
    }
    
    # Add on each sample to another column of tally.df
    tally.df[, each.sample] <- NA
    anno.idx <- interval_overlap(Intervals(sample.int[[each.sample]][,c(1,2)]), ref.interval)
    for(each.idx in c(1:length(anno.idx))){
      if(ncol(sample.int[[each.sample]]) == 3){
        tally.df[unlist(anno.idx[each.idx]), each.sample] <- sample.int[[each.sample]][each.idx,3]
      }
    }
  }
  return(tally.df)
}

# Uses rle() to obtain consecutive values in a dataframe and index them
getRleIdx <- function(x, col.id, na.val=-100){
  
  if(length(col.id) > 1) {
    uniq.id <- apply(x, 1, function(y) paste(y[col.id], collapse="-"))
    x$uniq <- uniq.id
    col.id <- 'uniq'
  }
  reformat.na.x <- as.character(x[,col.id])
  reformat.na.x[which(is.na(reformat.na.x))] <- na.val
  rle.x <- rle(reformat.na.x)

  
  #Get the array index for the start-to-end of each unique value/changepoint
  rle.x$start.idx <- c(1, (cumsum(rle.x$lengths) + 1)[-length(rle.x$lengths)])
  rle.x$end.idx <- rle.x$start.idx + (rle.x$lengths - 1)
  rle.x$values[which(rle.x$values == na.val)] <- NA
  rle.x$na.stat <- !is.na(rle.x$values)
  return(rle.x)
}

# Using the RLE outlined from getRleIdx, it collapses the segments and removes
collapseSegments <- function(x, cn.stat=TRUE, group.collapse=TRUE, group.ids=NA){
  # If cn.stat==FALSE, collapse based on LOH status rather than the Copy-number state
  if(!cn.stat){
    col.ids <- colnames(x)[-c(1:7)]
    x.cn.raw <- x
    x[,col.ids][!is.na(x[,col.ids])] <- 1
  }
  
  # If group.collapse==TRUE, collapsed based only a subset of samples
  if(group.collapse & !all(is.na(group.ids))){
    group.idx <- which(colnames(x) %in% group.ids)
    x.group.raw <- x
    x <- x[,c(1:7, group.idx)]
    x <- x[,c(1:3,5)]
  }
  
  # Get index of contiguous segments that need to be collapsed
  all.cols <- colnames(x)
  all.cols <- all.cols[which(!all.cols %in% c('start.pos', 'end.pos', 'length'))]
  x.rle <- getRleIdx(x, all.cols)
  
  # Restore original data-frame to be collapsed
  if(group.collapse & !all(is.na(group.ids))) x <- x.group.raw
  if(!cn.stat) x <- x.cn.raw
  
  x.list <- lapply(c(1:length(x.rle$values)), function(each.value){
    #Identify the start and end-cytoband for the given "ab" copy-state
    x.row <- x[x.rle$start.idx[each.value], , drop=FALSE]
    x.row[,'end.pos'] <- x[x.rle$end.idx[each.value], 'end.pos']
    x.row[,'length'] <- x.row[,'end.pos'] - x.row[,'start.pos']
    return(x.row)
  })
  x.df <- do.call("rbind", x.list)
  #x.df <- x.df[-which(x.df$length == 1),]
  if(any(colnames(x.df) == 'uniq')) x.df <- x.df[,-which(colnames(x.df) %in% 'uniq')]
  return(x.df)
}


compactLowRepeatSeg <- function(x){
  colnames(x)[c(1:3)] <- c("chr", "start.pos", "end.pos")
  x.split <- split(x, x$chr)
  compact.all.x <- lapply(names(x.split), function(each.chr){
    print(each.chr)
    if(nrow(x.split[[each.chr]]) > 0){
      x.interval <- Intervals(x.split[[each.chr]][, c('start.pos', 'end.pos')])
      
      x.union <- interval_union(x.interval)
      interval_overlap(x.interval, x.union)
      x.overlap.union <- interval_overlap(x.union, x.interval)
      
      compact.x.list <- lapply(c(1:length(x.overlap.union)), function(each.overlap){
        x.range <- x.overlap.union[[each.overlap]]
        x.start.row <- x.split[[each.chr]][x.range[1],, drop=FALSE]
        x.end.row <- x.split[[each.chr]][x.range[length(x.range)],,drop=FALSE]
        x.start.row$end.pos <- x.end.row$end.pos
        return(x.start.row)
      })
      compact.x <- do.call("rbind", compact.x.list)
      return(compact.x)
    }
  })
  x <- as.data.frame(do.call("rbind", compact.all.x), stringsAsFactors=FALSE)
  return(x)
}

