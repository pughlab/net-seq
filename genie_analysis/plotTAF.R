####################
######  Functions
{
  getPoly <- function(low.lim, up.lim, af, type){
    if(type %in% 'x'){
      poly.coord <- c(low.lim, af, up.lim, af, low.lim)
    } else if (type %in% 'y'){
      poly.coord <- c(0, 0.5, 0, -0.5, 0)
    }
    return(poly.coord)
  }
}

library(scales)
source("~/git/net-seq/genie_analysis/src/afCalculation.R")

####################
######  Setup
{
  pch.key <- c("ATRX" = 22,
               "DAXX" = 22,
               "MEN1" = 23,
               "Other"=16)
  col.key <- c("ATRX" = 'red',
               "DAXX" = 'red',
               "MEN1" = 'red',
               'Other'=alpha('black', 0.5))
  model.pos.key <- c("MEN1"=-1,
                     "DAXX"=-0.75,
                     "ATRX"=-0.5)
  purity.x <- -0.25
  cex.val <- 0.6
  rm.net.samples <- paste0("GENIE_Pnet_", c(16, 18, 19))
  rm.net.samples <- paste0("GENIE_Pnet_", c(23, 24, 26, 27, 30, 35))
  
  id.mapping <- read.table("~/git/net-seq/genie_analysis/data/genie_mapping.txt", 
                           sep="\t", header=FALSE)
  id.mapping <- lapply(split(id.mapping, id.mapping$V1), function(x) as.character(x[[2]]))
}

####################
######  Pre-process
{
  data.df <- read.table("~/git/net-seq/genie_analysis/data/Genie_purity_estimations.TAF2.txt",
                        sep="\t", header=TRUE, stringsAsFactors = FALSE, check.names = FALSE)
  
  # Filter data from the dataset
  data.df <- data.df[which(data.df$purity > 0.25),]
  data.df <- data.df[which(data.df$gene != 'ATRX-Ploidy'),]
  
  # EStablish the plotting order
  otb.ord <- c('29', '05', '30', '31', '34', '25', '03', '13', '17', '02', '11', '15', '32')
  net.ord <- c(paste0("NET-", c("001t2", "003t1", "003t2", "008t2", "009t1", "009t2")),
               paste0("GENIE_Pnet_", c(9:37)),
               paste0("NET-1", otb.ord))
  net.ord <- net.ord[-which(net.ord %in% rm.net.samples)]
  
  # Generate sample-specific data
  data.df$uid <- paste0(data.df$sample.id, data.df$extra.tag)
  data.split <- split(data.df, data.df$uid)
  data.split <- data.split[rev(net.ord)]
  data.split <- data.split[-which(is.na(names(data.split)))]
  
  
  taf.plot.vals <- lapply(data.split, function(sample.df, pch.key, model.pos.key){
    # Get the theoretical allelic fraction estimates (Ft)
    taf.vals <- apply(sample.df, 1, function(x) calcTAF(Falt=x['obs.af'], Pn=x['Pn'], pur.t=x['purity'], Pt=x['Pt']))
    falt.vals <- apply(sample.df, 1, function(x) calcFalt(Ft=1, Fn=0, Pn=x['Pn'], pur.t=x['purity'], Pt=x['Pt']))
    pur.vals <- apply(sample.df, 1, function(x)  calcPurity(Ft=1, Fn=0, Pn=x['Pn'], Falt=x['obs.af'], Pt=x['Pt']))
    pur.vals <- rep(mean(pur.vals), length(pur.vals))
    
    # Get the y.position based on sample id index
    sample.id <- unique(sample.df$uid)
    sample.val <- rep(match(sample.id, names(data.split)), length(taf.vals))
    
    # Get the right PCH keys for ATRX/DAXX and MEN1
    gene.ids <- sample.df$gene
    gene.ids.idx <- unlist(sapply(c("MEN1", "DAXX", "ATRX"), function(x) grep(x, gene.ids)))
    gene.ids[-gene.ids.idx] <- 'Other'
    gene.ids <- gsub("_snp", "", gene.ids)
    pch.vals <- pch.key[gene.ids]
    col.vals <- col.key[gene.ids]
    
    # Gets the ploidy and purity models
    model.vals <- paste(sample.df$Balt, sample.df$Pt, sep="/")
    model.text.vals <- model.pos.key[sample.df$gene]
    purity.vals <- rep(mean(as.numeric(as.character(sample.df$purity))), nrow(sample.df))

    
    return(data.frame("y.pos"=sample.val,
                      "x.pos"=taf.vals,
                      "pch.vals"=pch.vals,
                      "col.vals"=col.vals,
                      "model.vals"=model.vals,
                      "model.text.vals"=model.text.vals,
                      "purity.vals"=purity.vals,
                      "sample.id"=sample.df$sample.id,
                      "gene.id"=sample.df$gene))
  }, pch.key=pch.key, model.pos.key=model.pos.key)
  
  taf.plot.vals <- do.call("rbind", taf.plot.vals)
  taf.plot.vals <- taf.plot.vals[which(!is.na(taf.plot.vals$x.pos)),]
}
  
####################
######  Main/Visualization
dir.create("~/Desktop/netseq/genie/plots/", recursive = TRUE, showWarnings = FALSE)
pdf("~/Desktop/netseq/genie/plots/theorAfEstimation_net-otb-genie.pdf", height=6.5, width = 7)
{
  par(mar=c(5.1, 9.1, 4.1, 2.1))
  with(taf.plot.vals, plot(x=x.pos, y=y.pos, pch=pch.vals, axes=FALSE, ylab='', xlab='',
                           xlim=c(-1, 1.7), col=as.character(col.vals)))
  abline(v = c(0.5, 1), lty=3, col="black")
  row.labels <- unique(gsub("\\..*", "", rownames(taf.plot.vals)))
  row.labels <- gsub("t1", "_T1", row.labels)
  row.labels <- gsub("t2", "_T2", row.labels)
  row.labels <- as.character(id.mapping[row.labels]) # Mapping
  axis(side = 2, at = unique(taf.plot.vals$y.pos), labels = row.labels,
       las=2, cex.axis=cex.val, tick=FALSE)
  # axis(side = 2, at = unique(taf.plot.vals$y.pos), labels = c(paste0("NET_1", seq(3, 1, by=-1)),
  #                                                             paste0("GENIE_", seq(24, 1, by=-1)),
  #                                                             paste0("NETseq_", seq(5, 1, by=-1))),
  #                                                             las=2, cex.axis=cex.val, tick=FALSE)
  axis(side = 1, at = c(0, 0.5, 1, 1.5), labels=c(0, 0.5, 1, 1.5))
  axis(side = 3, at = c(model.pos.key, purity.x), labels=c(names(model.pos.key), "purity"), 
       cex.axis=cex.val, tick=FALSE, font=2)
  with(taf.plot.vals, text(x = model.text.vals, y=y.pos, labels = model.vals, adj=0, cex=cex.val))
  
  with(taf.plot.vals, segments(x0 = rep(purity.x - 0.1, nrow(taf.plot.vals)), 
                               y0 = y.pos, 
                               x1 = rep(purity.x + 0.1, nrow(taf.plot.vals)), 
                               y1= y.pos))
  conv.pur <- (taf.plot.vals$purity.vals / 5 ) - 0.1
  with(taf.plot.vals, points(x=purity.x + conv.pur, y=y.pos, col='black', pch=16))
  # with(taf.plot.vals, text(x = rep(purity.x, nrow(taf.plot.vals)), y=y.pos, 
  #                          labels = format(as.numeric(as.character(purity.vals)), 2), adj=0, cex=cex.val))
} 
dev.off()

pdf("~/Desktop/netseq/genie/plots/theorAfEstimation_legend.pdf", height=3, width = 7)
plot(0, type='n', xlim=c(1,10), ylim=c(0,6), xaxt='n', yaxt='n', ylab='', xlab='', axes=FALSE)
points(x=rep(1, 3), y=seq(1,5, by=2), pch=c(16, 22, 23), col=c(alpha('black', 0.5), "red", "red"))
text(x = rep(1.5, 3), y=seq(1,5, by=2), adj = 0, labels = c("Other", "ATRX/DAXX", "MEN1"), cex=0.7)
dev.off()


# Fishers exact test for enrichment
{
  hom.lim <- 0.85 # 0.63
  gene.data.list <- split(taf.plot.vals, taf.plot.vals$gene.id)
  
  binom.pval <- lapply(gene.data.list, function(gene.data.df) {
    num.of.samples <- nrow(gene.data.df)
    gene.data.df$hom.stat <- gene.data.df$x.pos >= hom.lim
    tryCatch({
      binom.val <- binom.test(x=sum(gene.data.df$hom.stat, na.rm=TRUE), 
               n=length(na.omit(gene.data.df$hom.stat)), 
               p = 0.5, alternative = c("greater"),
               conf.level = 0.95)$p.val
      c(binom.val, length(na.omit(gene.data.df$hom.stat)))
    }, error=function(e) { c(NA, NA) })
  })
  binom.pval <- do.call("rbind", binom.pval[c("MEN1", "DAXX", "ATRX")])
  p.adjust(binom.pval[,1], method = "bonferroni")
  
  
}
