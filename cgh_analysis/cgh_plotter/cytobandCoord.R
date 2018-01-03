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



#### Get's the cytoband locations for a set of genomic locations
getCytobands <- function(net.sequenza.df, cytoband.df, 
                         ret.stat='filter', collapseAB=FALSE){
  # If collapseAB is TRUE, it will collapse all gains regardless of degree of copy-state
  if(collapseAB){
    gain.idx <- which(net.sequenza.df$CNt > 2)
    gain.loh.idx <- which(with(net.sequenza.df[gain.idx,], A == 0 | B == 0))
    
    net.sequenza.df[gain.idx,]$A <- 3
    net.sequenza.df[gain.idx,]$B <- 1
    if(length(gain.loh.idx) > 0) net.sequenza.df[gain.loh.idx,]$B <- 0
  }
  
  # Creates unique rows of copy-states and splits based on chromosome
  net.sequenza.df$ab <- paste(net.sequenza.df$A, net.sequenza.df$B, sep="-")
  net.sequenza.list <- split(net.sequenza.df, net.sequenza.df$chromosome)
  cytoband.list <- split(cytoband.df, cytoband.df$chrom)
  
  
  cn.mdat.chr <- list()
  cn.mdat.filt.chr <- list()
  for(each.chr in paste("chr", c(1:22, "X", "Y"), sep="")){
    # Identifies overlap between CN-segments and cytoband intervals
    net.int <- Intervals(net.sequenza.list[[each.chr]][,c('start.pos', 'end.pos'), drop=FALSE])
    cytoband.int <- Intervals(cytoband.list[[each.chr]][,c('chromStart', 'chromEnd'), drop=FALSE])
    int.overlap <- interval_overlap(net.int, cytoband.int)
    length.intersection <- interval_intersection(net.int, cytoband.int)
    
    cn.mdat.list <- lapply(c(1:nrow(net.int)), function(each.interval){
      #overlapping segments from cytoband and CN-segment
      cytoband.mdat <- cytoband.list[[each.chr]][int.overlap[[each.interval]],,drop=FALSE]
      cn.mdat <- net.sequenza.list[[each.chr]][each.interval,,drop=FALSE]
      
      if(nrow(cytoband.mdat) == 1){
        # If overlap with only 1 cytoband interval, calculate size of overlap
        length.val <- length.intersection[each.interval,2] - length.intersection[each.interval,1]
        band.size <- with(cytoband.mdat, chromEnd - chromStart)
        band.fragment <- round(length.val / band.size, 2)
        band.name <- cytoband.mdat$name
      } else {
        # If overlap with multiple cytoband interval, concatenate the range
        band.fragment <- 1
        band.name <- cytoband.mdat$name[1]
        band.name <- paste(band.name, 
                           cytoband.mdat$name[nrow(cytoband.mdat)], sep="-")
      }
      add.df <- data.frame("cytoband"=band.name,
                           "band.fraction"=band.fragment)
      return(cbind(cn.mdat, add.df))
    })
    ## > Creates the unfiltered cytoband labelled segments
    cn.mdat.chr[[each.chr]] <- do.call("rbind", cn.mdat.list)
    
    # Collapse cytoband location for regions of the genome that don't change in copy-state
    cn.mdat.filt <- cn.mdat.chr[[each.chr]][which(cn.mdat.chr[[each.chr]]$band.fraction > 0.25),] # Remove small artifactual segments
    cytoband.rle <- getRleIdx(cn.mdat.filt, 'ab')
    cn.mdat.cytoband <- lapply(c(1:length(cytoband.rle$values)), function(each.value){
      #Identify the start and end-cytoband for the given "ab" copy-state
      cn.mdat.row <- cn.mdat.filt[cytoband.rle$start.idx[each.value], , drop=FALSE]
      cn.mdat.row[,'end.pos'] <- cn.mdat.filt[cytoband.rle$end.idx[each.value], 'end.pos']
      start.cytoband <- cn.mdat.filt[cytoband.rle$start.idx[each.value],'cytoband']
      start.cytoband <- strsplit(as.character(start.cytoband), split="-")[[1]][1]
      end.cytoband <- cn.mdat.filt[cytoband.rle$end.idx[each.value],'cytoband']
      end.cytoband <- strsplit(as.character(end.cytoband), split="-")[[1]]
      end.cytoband <- end.cytoband[length(end.cytoband)]
      
      # Format the reporting cytoband value appropriately
      cytoband.val <- paste(start.cytoband, end.cytoband, sep="-") 
      if(start.cytoband == end.cytoband) cytoband.val <- start.cytoband
      
      cn.mdat.row[,'cytoband'] <- cytoband.val
      return(cn.mdat.row)
    })
    ## > Creates the fitlered cytoband labelled segments
    cn.mdat.filt.chr[[each.chr]] <- do.call("rbind", cn.mdat.cytoband)
  }
  cn.mdat <- do.call("rbind", cn.mdat.chr)
  cn.mdat.filt <- do.call("rbind", cn.mdat.filt.chr)
  
  
  
  if(ret.stat == 'filter'){
    return(cn.mdat.filt)
  } else if(ret.start == 'unfilter'){
    return(cn.mdat)
  } else {
    print("Please indicate ret.start as either 'filter' or 'unfilter'")
  }
}

# Summarizes the dataframe into cytoband loci for neutral, gains and losses (LOH/non-LOH)
getCytobandSummary <- function(net.cytoband.df){
  summ.cytoband <- list()
  summ.cytoband[['neutral']] <- summarizeCytoband(net.cytoband.df[which(net.cytoband.df$CNt == 2),])
  summ.cytoband[['gain']] <- summarizeCytoband(net.cytoband.df[which(net.cytoband.df$CNt > 2),])
  summ.cytoband[['loss']] <- summarizeCytoband(net.cytoband.df[which(net.cytoband.df$CNt < 2),])
  
  summ.cytoband.df <- do.call("rbind", summ.cytoband)
  return(summ.cytoband.df)
}


# Splits the cytoband information into LOH and Non-LOH to collapse down into one liners
summarizeCytoband <- function(subset.cyto.df){
  loh.idx <- with(subset.cyto.df, A != 0 & B != 0)
  loh.df <- subset.cyto.df[which(!loh.idx),]
  nonLoh.df <- subset.cyto.df[which(loh.idx),]
  
  collapseCytoband <- function(x){
    if(nrow(x) == 1){
      unique(x$chromosome)
    } else {
      paste(unique(x$chromosome), x$cytoband, sep='')
    }
  }
  cytoband.het.list <- lapply(split(nonLoh.df, nonLoh.df$chromosome), collapseCytoband)
  cytoband.loh.list <- lapply(split(loh.df, loh.df$chromosome), collapseCytoband)
  
  return(list('loh'=unlist(cytoband.loh.list),
              'het'=unlist(cytoband.het.list)))
}
