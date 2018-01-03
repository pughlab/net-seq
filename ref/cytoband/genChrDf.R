# Uses rle() to obtain consecutive values in a dataframe and index them
getRleIdx <- function(x, col.id, na.val=-100){
  reformat.na.x <- x[,col.id]
  reformat.na.x[which(is.na(reformat.na.x))] <- na.val
  rle.x <- rle(reformat.na.x)
  
  #Get the array index for the start-to-end of each unique value/changepoint
  rle.x$start.idx <- c(1, (cumsum(rle.x$lengths) + 1)[-length(rle.x$lengths)])
  rle.x$end.idx <- rle.x$start.idx + (rle.x$lengths - 1)
  rle.x$values[which(rle.x$values == na.val)] <- NA
  rle.x$na.stat <- !is.na(rle.x$values)
  return(rle.x)
}

cytoband.df <- read.csv("/Users/rquevedo/git/reference/cytoband/cytoband.hg19.tsv", 
                        header = TRUE, sep="\t", stringsAsFactors = FALSE, check.names = FALSE)

colnames(cytoband.df)[1] <- 'chrom'
chr.rle <- getRleIdx(cytoband.df, 'chrom')
chr.order <- paste("chr", c(1:22, "X", "Y"), sep="")

# Subset the cytoband.df for just Chr 1-22,X,Y start/end coordinates and generate cumulative lengths
cum.start <- 0; start.end.length <- 0
chr.df <- data.frame()
for(each.chr in chr.order){
  rle.idx <- match(each.chr, chr.rle$values)
  
  cum.start <- cum.start + start.end.length 
  start.end.length <- with(cytoband.df, chromEnd[chr.rle$end.idx[rle.idx]], 
                           chromStart[chr.rle$start.idx[rle.idx]])
  cum.end <- cum.start + start.end.length
  
  temp.chr.df <- data.frame("chrom" = cytoband.df[chr.rle$start.idx[rle.idx],'chrom'],
                            'chromStart' = cytoband.df[chr.rle$start.idx[rle.idx],'chromStart'],
                            'chromEnd' = cytoband.df[chr.rle$end.idx[rle.idx],'chromEnd'],
                            'cumStart' = cum.start,
                            'cumEnd' = cum.end)
  chr.df <- rbind(chr.df, temp.chr.df)
}
chr.df$chr <- gsub("chr", "", chr.df$chrom)

# Generate the same cumulative start-end Dataframe for cytoband level dataframe of just chr1-22, X, Y
cum.start <- 0; start.end.length <- 0
raw.cytoband.df <- cytoband.df
cytoband.list <- split(cytoband.df, cytoband.df$chrom)
for(each.chr in chr.order){
  cytoband.list[[each.chr]]$length <-  with(cytoband.list[[each.chr]], chromEnd - chromStart)
  cytoband.list[[each.chr]]$cumStart <- cum.start + c(0, cumsum(cytoband.list[[each.chr]]$length)[-nrow(cytoband.list[[each.chr]])])
  cytoband.list[[each.chr]]$cumEnd <- cum.start + cumsum(cytoband.list[[each.chr]]$length)
  cum.start <- cytoband.list[[each.chr]][nrow(cytoband.list[[each.chr]]),'cumEnd']
}
cytoband.df <- do.call("rbind", cytoband.list[chr.order])

save(cytoband.df, raw.cytoband.df, chr.df, 
     file="/Users/rquevedo/git/reference/cytoband/chr_cytoband.hg19.Rdata")
