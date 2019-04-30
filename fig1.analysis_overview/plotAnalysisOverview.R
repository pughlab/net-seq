## Creates a circle-line plot for analysis done on NETseq
#### Parameters ####
line.width <- 0.055 # Size of the line going across platforms
point.size <- 2 # Size of the dots for identity matrix
text.size <- 0.7
dataset.idx <- 4  # y loci for the horizontal line dividing dataset from platforms

#### Functions ####
blankPlot <- function(xlim, ylim){
  plot(0, type='n', xlab='', ylab='',
       xaxt='n', yaxt='n', axes=FALSE, 
       xlim=xlim, ylim=ylim)
}

#### MAIN ####
ov <- read.table("~/git/net-seq/fig1.analysis_overview/overview.csv",
                 sep=",", header=TRUE, stringsAsFactors = FALSE)
ov <- ov[,-which(apply(ov, 2, function(x) all(is.na(x))))]
ov.IM <- t(ov[,-c(1:3)])
colnames(ov.IM) <- ov$platform
ov.IM <- ov.IM[nrow(ov.IM):1,]




pdf("~/git/net-seq/fig1.analysis_overview/overview.pdf", height = 7, width=5)
split.screen(matrix(c(0, 0.7, 0.6, 1.0,
                      0, 0.7, 0, 0.6,
                      0.7, 1, 0, 0.6), byrow = TRUE, ncol=4))
### Plot identity matrix
screen(2); par(mar=c(1, 2, 0.5, 0.2));
n.anal <- nrow(ov.IM)
n.plat <- ncol(ov.IM)
blankPlot(xlim=c(0, n.plat), ylim=c(1, n.anal))

sapply(1:n.anal, function(i) {
  ## Background tracks
  points(x=c(1:n.plat), y=rep(i, n.plat), 
         pch=16, col="grey", cex=point.size)
  
  ## Foreground tracks based on identity matrix
  rect.track <- which(as.logical(ov.IM[i,]))
  rect(xleft = min(rect.track), ybottom = i-line.width, 
       xright = max(rect.track), ytop = i+line.width,
       lty =, col="black", border=FALSE)
  points(x=rect.track, y=rep(i, length(rect.track)),
         pch=16, col="black", cex=point.size)
})


### Label analysis methods
screen(3); par(mar=c(1, 0, 0.5, 0));
blankPlot(xlim=c(0, 10), ylim=c(1, n.anal))
anal.idx <- seq_along(rownames(ov.IM))
text(y=anal.idx, x=0, pos = 4,
     labels = gsub("\\.", " ", rownames(ov.IM)),
     cex=text.size)
abline(v = 0, lty=1, lwd=2, col="black")

### Label platforms used
screen(1); par(mar=c(0, 2, 0, 0.2));
blankPlot(xlim=c(0, n.plat), ylim=c(0, 10))
plat.idx <- seq_along(colnames(ov.IM))
text(x=plat.idx, y=0, adj=0, srt=90,
     labels = colnames(ov.IM),
     cex=text.size)
text(x=plat.idx, y=dataset.idx, adj=0, srt=90,
     labels = ov$dataset,
     cex=text.size)

for(i in ov$dataset[ov$dataset != '']){
  s.idx <- match(i, ov$dataset) - 0.4

  seg.len <- ov$dataset[match(i, ov$dataset):nrow(ov)] == ''
  seg.len[1] <- TRUE
  e.idx <- (rle(seg.len)$lengths[1] + match(i, ov$dataset)) - 0.6
  
  rect(xleft = s.idx, ybottom = (dataset.idx - 0.8), 
       xright = e.idx, ytop = (dataset.idx - 0.6),
       col="black", border=FALSE)
}


close.screen(all.screens=TRUE)
dev.off()