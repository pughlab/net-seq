#GTEx Visualization Functions

# Plots the coloured boxplots used in gtex.geneTpmPlotter
plotTissueBoxplots <- function(tissue.list, min.val=-1, max.val, tissue.col, col.n='log2(TPM + 1)'){
  boxplot(tissue.list, ylim=c(min.val, max.val), 
          ylab=col.n, las=2, cex.axis=0.5)
  tissue.cnt <- 1
  for(each.tissue in tissue.list){
    points(rep(tissue.cnt, length(each.tissue)), 
           each.tissue, 
           col=alpha(tissue.col[tissue.cnt], 0.4), pch=16)
    tissue.cnt <- tissue.cnt + 1
  }
  boxplot(tissue.list, ylim=c(min.val, max.val), xaxt='n', yaxt='n', add=TRUE)
}

# Plots the coloured 95% CI rectangles for Tissue-specific TPM
plot95Ci <- function(tissue.list, min.val=0, max.val=100, tissue.col, goi=NULL){
  # Set up blank plot
  plot(0, type="n", xlim=c(0, length(tissue.list)), ylim=c(min.val,max.val), 
       axes=FALSE, ylab="TPM", xlab="",  main=goi)
  axis(side = 2, at = seq(min.val, max.val, by=((max.val - min.val)/3)), 
       labels = seq(min.val, max.val, by=((max.val - min.val)/3)), las=2,
       cex.axis=0.5, cex=0.8)
  
  # Add in confidence interval rectangles
  for(each.goi.idx in c(1:length(tissue.list))){
    ci.vals <- tissue.list[[each.goi.idx]]
    tissue.id <- gsub(" - .*$", "", names(tissue.list)[each.goi.idx])
    rect(xleft=(each.goi.idx - 0.2), ybottom=ci.vals[2],
         xright=(each.goi.idx + 0.2), ytop=ci.vals[3], border=FALSE, 
         col=tissue.col.ord[tissue.id])
    rect(xleft=(each.goi.idx - 0.2), ybottom=(ci.vals[1] - 0.5),
         xright=(each.goi.idx + 0.2), ytop=(ci.vals[1] + 0.5), border=FALSE, col="black")
  }
}


# Plots a tissue type associated with a given colour
plotTissueLegend <- function(tissue.col){
  plot(1, type='n', axes=FALSE, ylab="", xlab="",
       ylim=c(0, 35), xlim=c(0,10))
  rect(xleft=rep(0, length(tissue.col)), 
       ybottom=seq(34, (34 - length(tissue.col)), by=-1),
       xright=rep(1, length(tissue.col)), 
       ytop=seq(35, (35-length(tissue.col)), by=-1), col=tissue.col)
  text(x=rep(1, length(tissue.col)), 
       y=seq(34.5, (34.5-length(tissue.col)), by=-1), 
       pos=4, labels=names(tissue.col), cex=0.5)
}