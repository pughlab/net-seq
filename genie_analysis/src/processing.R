# Plots the variants on top of a CN-segment plot
plotPoints <- function(mut.x, cn.x, mut.col, snp.stat='dbsnp', increase.interval.size=FALSE){
  if(increase.interval.size){
    cn.x$loc.start <- cn.x$loc.start - 10
    cn.x$loc.end <- cn.x$loc.end + 10
  }
  cn.int <- Intervals(as.matrix(cn.x[,c("loc.start", "loc.end")]))
  mut.int <- as.matrix(mut.x[,c("start_position", "end_position")])

  row.match.idx <- apply(mut.int, 1, function(mut.single){
    mut.single <- Intervals(mut.single)
    which(sapply(interval_included(cn.int, mut.single), function(x) length(x) > 0))[1]
  })
  
  if(snp.stat == 'novel') snp.stat <- 1 else snp.stat <- -1
  with(mut.x, points(x=start_position, y=cn.x[row.match.idx, 'seg.mean'], pch=16, col=mut.col, bg=mut.col))
  with(mut.x, text(x=start_position, 
                   y=(cn.x[row.match.idx, 'seg.mean'] - (snp.stat*seq(0.1, 0.1*length(row.match.idx), by=0.1))), 
                   labels=round(alt.af, 2), cex=0.7, col=mut.col))
  with(mut.x, text(x=with(cn.x, max(loc.end)/2), 
                   y=rev(snp.stat*(-1.5 + seq(0.1, 0.1*length(row.match.idx), by=0.1))), 
                   labels=hugo_symbol, cex=0.7, col=mut.col))
}



# Gets the predominant CN state for each set of CN-equivalent genes
getModeCnState <- function(x){
  # x is a dataframe containing all theorized purities and cn-states for all genes
  all.cn.vals <- seq(1:8)
  x.split <- split(x, x$gene, drop = FALSE)
  
  # Matches each genes entry to all possible copy-states
  x.match <- lapply(x.split, function(i, j){
    match(i$cn, j)
  }, j=all.cn.vals)
  
  # Identifies the copy-state that is conserved across most/all of the genes
  x.table <- table(unlist(x.match))
  cn.state.mode <- x.table[which(x.table == max(x.table))]
  return(cn.state.mode)
}