# Function: populateGenomeBins
# Purpose:  Takes in a pre-determined blank list of all bin sizes for each chromosome in the human genome
#   and maps out the given segment length and CN-segments to each of the bins.  This will allow for an easier
#   bin-vs-bin comparison in future steps.
# Input:  cn.val.f1: dataframe containing the CN segments and their seg-values (hscr.a1, seg_mean, etc..)
#         chr.bins: The blank list containing all preformatted empty bins
#         seg.col.name: column name for the segment header (i.e. seg_mean, hscr.a1, hscr.a2, etc...)
# Returns:  genome.tbins.df: Dataframe containing all the preformatted bins with their associated seg-values 
#   according to cn.val.f1.
populateGenomeBins <- function(cn.val.f1, chr.bins, seg.col.name){
  mismatch.val <- 0
  
  genome.tbins.df <- data.frame(matrix(ncol=4, nrow=0))
  for(each.chr in global.chrom){
    composite.bin.df <- as.data.frame(matrix(ncol=5, nrow=0),
                                      stringsAsFactors=FALSE, check.names=FALSE) #Keeps track of segments that do not fit whole bins
    matching.intervals.df <- as.data.frame(matrix(ncol=3, nrow=0),
                                           stringsAsFactors=FALSE, check.names=FALSE) # Keeps track of whole bin segments.
    
    colnames(matching.intervals.df) <- c("Start.bin", "End.bin", seg.col.name)
    blank.chr.int <- Intervals(chr.bins[[each.chr]], closed=TRUE)
    # Very sloppy QC checks for the input of CN-abberations and Start-end bp locations
    if(length(cn.val.f1) > 0){
      cn.chr.f1 <- cn.val.f1[cn.val.f1$Chromosome %in% each.chr,]
      if(dim(cn.chr.f1)[1] > 0){
        if(cn.chr.f1$End.bp == 0){
          cn.chr.f1 <- as.data.frame(matrix(ncol=5, nrow=0))
          colnames(cn.chr.f1) <- c("chr.chr", "Start.bp", "End.bp", seg.col.name, "Chromosome")
        }
      } 
    } else {
      cn.chr.f1 <- as.data.frame(matrix(ncol=5, nrow=0))
      colnames(cn.chr.f1) <- c("chr.chr", "Start.bp", "End.bp", seg.col.name, "Chromosome")
    }
    if(dim(cn.chr.f1)[1] > 0){          ### Checks if there is at least one CN instance in the current chromosome
      for(each.chr.row in 1:dim(cn.chr.f1)[1]){
        cn.int.f1 <- Intervals(as.numeric(matrix(c(as.integer(cn.chr.f1[each.chr.row, 'Start.bp']), 
                                                   as.integer(cn.chr.f1[each.chr.row, 'End.bp'])), ncol=2)))
        match.matrix <- as.matrix(interval_intersection(cn.int.f1, blank.chr.int))
        # Matches the blank matrix to the "match matrix" that has removed first and last row (to be calculated later)
        temp.matrix <- matrix(sort(intersect(chr.bins[[each.chr]], match.matrix[-c(1,nrow(match.matrix)),])), ncol=2, byrow=T)
        temp.matrix <- cbind(temp.matrix, rep(cn.chr.f1[each.chr.row, seg.col.name], nrow(temp.matrix)))
        colnames(temp.matrix) <- c("Start.bin", "End.bin", seg.col.name)
        matching.intervals.df <- rbind(matching.intervals.df, as.data.frame(temp.matrix, stringsAsFactors=FALSE))
        colnames(matching.intervals.df) <- c("Start.bin", "End.bin", seg.col.name)
        
        # Stores the intervals that do not fill a complete bin in the composite dataframe to handle after.
        # bp-level intervals match in the bin
        first.match.row <- match.matrix[1,]
        last.match.row <- match.matrix[nrow(match.matrix),]
        # the bin being matched to
        first.match.blank.int <- Intervals(first.match.row)
        last.match.blank.int <- Intervals(last.match.row)
        composite.matrix <- matrix(c(first.match.row, blank.chr.int[!is.na(interval_overlap(blank.chr.int, first.match.blank.int) > 0)],
                                     last.match.row, blank.chr.int[!is.na(interval_overlap(blank.chr.int, last.match.blank.int) > 0)]), 
                                   nrow=2, byrow=T)
        composite.matrix <- cbind(composite.matrix, rep(cn.chr.f1[each.chr.row, seg.col.name], nrow(composite.matrix)))
        if(length(which(first.match.row == last.match.row)) ==  2){   # If the first.match.row vector matches last.match.row vector
          composite.bin.df <- data.frame(rbind(as.matrix(composite.bin.df), composite.matrix[1,]),
                                         stringsAsFactors=FALSE, check.names=FALSE)
        } else {
          composite.bin.df <- data.frame(rbind(as.matrix(composite.bin.df), composite.matrix),
                                         stringsAsFactors=FALSE, check.names=FALSE)
        }
        # Composite matrix result: start.bp, end.bp, start.bin, end.bin, segment.col.name
        
      }
      colnames(composite.bin.df) <- c("Start.bp", "End.bp", "Start.bin", "End.bin", seg.col.name)
      matching.intervals.df <- rbind(matching.intervals.df, generateWeightedCompositeBin(composite.bin.df, seg.col.name))
      matching.intervals.df <- matching.intervals.df[order(matching.intervals.df$Start.bin),]
      #Finds the intervals not included in the analysis and fills them in
      match.intervals <- Intervals(as.numeric(matrix(c(matching.intervals.df[, 'Start.bin'], 
                                                       matching.intervals.df[, 'End.bin']), ncol=2, byrow=F)))
      mismatched.int.matrix <- as.matrix(interval_difference(blank.chr.int, match.intervals))
      mismatched.int.matrix <- cbind(mismatched.int.matrix, rep(mismatch.val, dim(mismatched.int.matrix)[1]))
      matching.intervals.df <- data.frame(rbind(as.matrix(matching.intervals.df), mismatched.int.matrix),
                                          stringsAsFactors=FALSE, row.names=NULL)
      
      
      matching.intervals.df <- matching.intervals.df[order(matching.intervals.df$Start.bin),]
      matching.intervals.df$chr <- each.chr
      
      
      
      genome.tbins.df <- rbind (genome.tbins.df, matching.intervals.df)
      # Need to find the complement of the intervals using interval_complement
    } else {        ### Fills in a mismatch.val (CN-state of 0) 
      matching.intervals.df <- as.data.frame(as.matrix(blank.chr.int))
      matching.intervals.df[,seg.col.name] <- 0
      matching.intervals.df$Chr <- each.chr
      colnames(matching.intervals.df) <- c("Start.bin", "End.bin", seg.col.name, "chr")
      
      genome.tbins.df <- rbind (genome.tbins.df, matching.intervals.df)
    }
    
  }
  colnames(genome.tbins.df) <- c("Start.bin", "End.bin", seg.col.name, "Chr")
  return(genome.tbins.df)
}

# Function: generateWeightedCompositeBin
# Purpose:  Looks at predetermined bins and proportionally weighs each CN-segment that makes up that bin
#   to return a weighted-CN-segment for that bin.
# Input:  composite.bin.df <- dataframe containing the actual intervals that make up each predetermined bins
#         seg.col.name <- column name for the segment header (i.e. seg_mean, hscr.a1, hscr.a2, etc...)
# Returns:  composite.bin.final.df: data-frame containing the genome binned into the set sizes and having a 
#   proportionally weighted CN segment according to to the segment length in that bin.
generateWeightedCompositeBin <- function(composite.bin.df, seg.col.name){
  composite.bin.final.df <- as.data.frame(matrix(ncol=3, nrow=0), 
                                          stringsAsFactors=FALSE, check.names=FALSE)
  composite.bin.df$Start.bin <- as.integer(composite.bin.df$Start.bin)
  composite.bin.df$End.bin <- as.integer(composite.bin.df$End.bin)
  for(each.uniq.bin in unique(composite.bin.df$Start.bin)){
    unique.bin.df <- composite.bin.df[composite.bin.df$Start.bin == each.uniq.bin, ]
    
    start.bin <- each.uniq.bin
    end.bin <- unique.bin.df[1,'End.bin']
    
    # Calculates the proportionally-weighted CN segment for the given bins    
    unique.bin.df$length <- as.integer(unique.bin.df$End.bp) - as.integer(unique.bin.df$Start.bp)
    total.cov <- sum(unique.bin.df$length)
    unique.bin.df$weighted.seg <- (unique.bin.df$length / total.cov) * as.numeric(unique.bin.df[,seg.col.name])
    
    composite.sum <- sum(unique.bin.df$weighted.seg)
    composite.bin.final.df <- rbind(composite.bin.final.df, 
                                    data.frame(start.bin, end.bin, composite.sum))
  }
  colnames(composite.bin.final.df) <- c("Start.bin", "End.bin", seg.col.name)
  
  return(composite.bin.final.df)
}

# Function: generateGenomeBins
# Purpose:  Uses the UCSC ChromInfo file to generate bins of a specific size for later on
#   CN evaluation
# Input:  - chr.df <- ucsc.chrom.info dataframe containing Chr and Size
#         - bin.s <- The size of the bins you want to create
# Returns:  genome.binned: data-frame containing the genome binned into the set sizes
generateGenomeBins <- function(chr.df, bin.s){
  # Formats the ucsc chrom info dataframe in an easier to use way.
  chr.df <- chr.df[ chr.df$chrom %in% global.chrom,]
  rownames(chr.df) <- chr.df$chrom
  chr.df <- chr.df[,-1]
  
  # Creates a List containing binned intervals for ezch chromsome
  chr.bin.list <- list()
  for(each.chr in global.chrom){
    bin.start.pos <- 1
    bin.max.matrix <- matrix(ncol=2, nrow=0)
    while(bin.start.pos < max(chr.df[as.character(each.chr),'size'])){
      bin.end.pos <- bin.start.pos + bin.s
      if(bin.end.pos > max(chr.df[as.character(each.chr),'size'])){
        bin.end.pos <- max(chr.df[as.character(each.chr),'size'])
      }
      bin.max.matrix <- rbind(bin.max.matrix, c(bin.start.pos, bin.end.pos))
      bin.start.pos <- bin.end.pos + 1
    }
    chr.bin.list[[as.character(each.chr)]] <- bin.max.matrix
  }
  
  return(chr.bin.list)
}




getSampleDf <- function(x, stat){
  sample.df <- as.data.frame(do.call("rbind", x), stringsAsFactors=FALSE)
  if(length(cnt.idx <- grep("CNt", colnames(sample.df))) > 0 ){
    sample.df <- sample.df[,-cnt.idx]
  }
  if(dim(sample.df)[2] > 0){
    # If the End position is smaller than the start position (e.g. p-arm of Chr 3)
    reord.row <- which(with(sample.df, as.integer(end) < as.integer(start)))
    if(length(reord.row) > 0){
      start.bp <- sample.df[reord.row, 'start']
      end.bp <- sample.df[reord.row, 'end']
      
      sample.df[reord.row, 'start'] <- end.bp
      sample.df[reord.row, 'end'] <- start.bp
    }
    
    sample.df$stat <- stat
    sample.df$chr <- gsub("^chr", "", sample.df$chr.chr)
    colnames(sample.df) <- c("chr.chr", "Start.bp", "End.bp", "stat", "Chromosome")
  }
  return(sample.df)
}



appendList <- function (x, val) 
{
  stopifnot(is.list(x), is.list(val))
  xnames <- names(x)
  if(!is.data.frame(x) && !is.data.frame(val)){
    for (v in names(val)) {
      x[[v]] <- if (v %in% xnames && is.list(x[[v]]) && is.list(val[[v]])){
        appendList(x[[v]], val[[v]])
      } 
    }
  } else {
    x <- rbind(x, val)
  }
  x
}


getCIFraction <- function(x, global.bin.size, chr.sel=c(1:22)){
  x <- x[which(x$Chr %in% chr.sel),]  # Filter out chr X, Y, MT
  #Generate bin sizes and identify ratio of bins without a stat=0 to total genome size analyed
  x$bin.size <- as.numeric(as.character(x$End.bin)) - as.numeric(as.character(x$Start.bin))
  x.ci <- x[which(x$stat != 0),]
  ci.frac <- sum(x.ci$bin.size) / sum(x$bin.size)
  return(ci.frac)
}

getCnSegs <- function(x.df, state){
  if(state %in% 'losses'){
    x.df <- tryCatch(x.df[which(x.df$CNt < 2),],
             error = function(e) {
               x.df <- NULL
             })
    id <- 'Losses'
    
  } else if(state %in% 'gains'){
    x.df <- tryCatch(x.df[which(x.df$CNt > 2),],
             error = function(e) {
               x.df <- NULL
             })
    id <- 'Gains'
  }else if(state %in% 'loh'){
    x.df <- x.df[which(x.df$A == 0 | x.df$B == 0),]
    id <- 'LOH'
  } else {
    warning("Improper state passed in")
    break
  }
  
  segs <- NULL
  if(nrow(x.df) > 0){
    rownames(x.df) <- paste(id, c(1:dim(x.df)[1]), sep="")
    segs <- apply(x.df, 1, function(x) list(c("chr.chr"=as.character(x['chromosome']), 
                                              "start"=as.character(x['start.pos']), 
                                              "end"=as.character(x['end.pos']),
                                              "CNt"=as.character(x['CNt']))))
    segs <- lapply(segs, function(x) unlist(x))
  } 
  return(segs)
}
