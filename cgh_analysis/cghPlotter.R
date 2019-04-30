library(gtools)
library(vcd)
library(intervals)
library(xtable)
require(RColorBrewer)

load("~/git/net-seq/ref/cytoband/chr_cytoband.hg19.Rdata")

hg19.cytoband <- read.csv("~/git/net-seq/ref/cytoband/cytoband.hg19.tsv", sep="\t", header=TRUE,
                          stringsAsFactors=FALSE, check.names=FALSE)
colnames(hg19.cytoband)[1] <- 'chrom'
load("~/git/loh_cn_visualize/data/data.aggregateLoh.Rdata")  #PNETs, GINETs, CCL
ucsc.chrom <- ucsc.chrom[order(ucsc.chrom$chr.num),]
ucsc.chrom$mid.st <- ucsc.chrom$chr.st.size + (ucsc.chrom$size / 2)
ucsc.chromInfo <- '~/Onedrive/bhk_lab/reference/ucsc.hg19.chromInfo.txt'
#load("~/git/loh_cn_visualize/data/data.aggregateLoh.Rdata")  #PNETs, GINETs, CCL
#load("~/git/shallowWgsCnv/data/otbData.aggregateLoh.medianLOH.Rdata")  # otb-PNET, otb-GINET
#load("~/git/shallowWgsCnv/data/otbData.aggregateLoh.muLOH.Rdata")  # otb-PNET, otb-GINET
load("~/Onedrive/PughLab/NET-seq/otb_samples/net_OTB/plots/postPathBlacklisted/medLOHseg-0.25/otbData.aggregateLoh.medianLOH.Rdata")  # Filtered one

source('~/git/loh_cn_visualize/cgh_plotter/binCgh.R')
source('~/git/loh_cn_visualize/cgh_plotter/cnaPlot.R')
source('~/git/loh_cn_visualize/cgh_plotter/genomicCoord.R')
source('~/git/loh_cn_visualize/cgh_plotter/cytobandCoord.R')
source('~/git/loh_cn_visualize/cgh_plotter/binLoh.R')
source("~/git/cn-annotater/src/annotateSeg.R")

cgh.dir <- '~/Onedrive/PughLab/NET-seq/cgh_analysis/all/'
setwd(cgh.dir)

eval.netseq <- 1  # Add in NETSEQ dataset

### ------------------------------ Pre-Processing ------------------------------ ###
### Converts the CGH data to genomic coordinates and pre-processes chromosome info
### ---------------------------------------------------------------------------- ###

### Reads in data and formats
cgh.colnames <- c("Study.id", "Case", "Sex", "Age", "NET", "F.stat", "Met", "Losses", "Gains")
cgh.list <- list()
for(each.file in list.files(pattern=".txt")){
  each.name <- gsub("\\..+$", "", each.file)
  cgh.list[[each.name]] <- read.table(each.file, sep="\t", header=TRUE,
                                      stringsAsFactors=FALSE, check.names=FALSE)
}
cgh.df <- do.call("smartbind", cgh.list)
cgh.df <- cgh.df[,cgh.colnames]
cgh.df[cgh.df$Losses %in% "None", 'Losses'] <- ""
cgh.df[cgh.df$Gains %in% "None", 'Gains'] <- ""


### Get genomic coordinates of CGH data and store in dataframe
gains.list <- apply(cgh.df, 1, function(x) getCytoPos(x['Gains'], hg19.cytoband))
loss.list <- apply(cgh.df, 1, function(x) getCytoPos(x['Losses'], hg19.cytoband))
cgh.df$Gains <-  gains.list
cgh.df$Losses <- loss.list

#Yes/No/Unknown Met status
cgh.df$Met.stat <- cgh.df$Met
cgh.df$Met.stat[which(!(cgh.df$Met %in% "No" | cgh.df$Met %in% "Unk"))] <- "Yes"

### Gets genomic information for plotting
t.genome.row <- which(ucsc.chrom$chr.num == 22)
t.genome.size <- ucsc.chrom[t.genome.row, 'chr.st.size']
t.genome.size <- t.genome.size + ucsc.chrom[t.genome.row, 'size']

### ------------------------------ Pre-Processing ------------------------------ ###
### Adding in NET-SEQ samples
### ---------------------------------------------------------------------------- ###
if(eval.netseq == 1){
  ### Gets all Losses and Gains from Sequenza
  net.seq.dir <- '~/Onedrive/PughLab/NET-seq/Sequenza/all'
  
  loss.list <- list()
  gain.list <- list()
  loh.list <- list()
  net.cytoband.list <- list()
  
  # populate a list for genomic regions of gains, loss and LOH based on 
  # Sequenza files for the discovery cohort
  for(each.net.file in list.files(path=net.seq.dir, pattern="T_2-N.gz_segments.txt")){
    each.name <- gsub("-T.+$", "", each.net.file)
    net.sequenza.df <- read.table(file.path(net.seq.dir, each.net.file), sep="\t", header=TRUE,
                                  stringsAsFactors=FALSE, check.names=FALSE)
    tryCatch({
    net.cytoband.df <- getCytobands(net.sequenza.df, cytoband.df, 
                                    ret.stat = 'filter', collapse=TRUE)
    net.cytoband.summary <- getCytobandSummary(net.cytoband.df)
    for(row.val in c(1:nrow(net.cytoband.summary))){
      for(col.val in c(1:ncol(net.cytoband.summary))){
        collapsed.val <- paste(unlist(net.cytoband.summary[row.val, col.val]), collapse=", ")
        net.cytoband.summary[row.val, col.val] <- collapsed.val
      }
    }
    net.cytoband.list[[each.net.file]] <- net.cytoband.summary
    
    loss.list[[each.name]] <- getCnSegs(net.sequenza.df, 'losses')
    gain.list[[each.name]] <- getCnSegs(net.sequenza.df, 'gains')
    loh.list[[each.name]] <- getCnSegs(net.sequenza.df, 'loh')},
    error=function(e){})
  }
  
  # populate a list for genomic regions of gains, loss and LOH based on 
  # Sequenza-like (made from sWGS tool) files for the validation cohort
  if(length(grep("otb", names(disease.list))) > 0){
    # names(disease.list)
    for(each.net in 'otb-pnet'){
      for(each.name in names(disease.list[[each.net]])){
        print(each.name)
        net.swgs.df <- as.data.frame(disease.list[[each.net]][[each.name]], stringsAsFactors=FALSE)
        colnames(net.swgs.df) <- c("chromosome", "start.pos", "end.pos", 
                                   "A", "B", "LOH")
        net.swgs.df$A <- as.integer(as.character(net.swgs.df$A))
        net.swgs.df$B <- as.integer(as.character(net.swgs.df$B))
        
        loh.idx <- which(net.swgs.df$LOH == 0)
        net.swgs.df[loh.idx,'B'] <- 1
        a.hi.idx <- intersect(which(net.swgs.df[,'A'] > 1), loh.idx)
        a.zero.idx <- intersect(which(net.swgs.df[,'A'] == 0), loh.idx)
        net.swgs.df[a.hi.idx,'A'] <- (net.swgs.df[a.hi.idx,'A'] - 1) 
        net.swgs.df[a.zero.idx,'A'] <- (net.swgs.df[a.zero.idx,'A'] + 1) 
        #net.swgs.df[which(net.swgs.df[loh.idx,'A'] > 1), 'A'] <- net.swgs.df[which(net.swgs.df[loh.idx,'A'] > 1), 'A'] - 1
        net.swgs.df$CNt <- with(net.swgs.df, as.integer(A) + as.integer(B))
        
        loss.list[[each.name]] <- getCnSegs(net.swgs.df, 'losses')
        gain.list[[each.name]] <- getCnSegs(net.swgs.df, 'gains')
        loh.list[[each.name]] <- getCnSegs(net.swgs.df, 'loh')
      }
    }
  }
  
  
  # Write out the cytoband coordinate files for NETseq.2017 samples
  dir.create(path = file.path(getwd(), "net_cytoband"), showWarnings = FALSE)
  for(each.name in names(net.cytoband.list)){
    write.table(net.cytoband.list[[each.name]],
                file=file.path("net_cytoband", 
                               paste(gsub(".gz.*", "", each.name), ".cytoband.txt", sep="")),
                quote=FALSE, col.names=TRUE, row.names=TRUE, sep="\t")
  }
  
  # Obtain the LOH overlap master list
  all.net.loh <- lapply(loh.list, function(x) as.data.frame(do.call("rbind", x))) #Combine all samples
  all.net.loh <- lapply(all.net.loh, function(x) split(x, x[,1])) #Split by Chromosome
  pnet.ids=c('NET-001', 'NET-003', 'NET-008', 'NET-009')
  ginet.ids=c('NET-002', 'NET-006', 'NET-007')
  populated.interval <- getLohOverlap(all.net.loh, 
                                      group.id1=c(names(disease.list[['otb-pnet']]), pnet.ids),
                                      group.id2=c(names(disease.list[['otb-ginet']]), ginet.ids))
  populated.interval <- collapseSegments(populated.interval,
                                         cn.stat=FALSE, group.collapse=TRUE,
                                         group.ids=c(names(disease.list[['otb-pnet']]), pnet.ids))
  source("~/git/cn-annotater/src/annotateSeg.R")
  populated.interval <- annotateSeg(populated.interval,
                                    loc.headers=c("chr", "start.pos", "end.pos"))
  
  populated.n1.interval <- populated.interval[which(populated.interval$`PNETs_(n=17)` <= 1),]
  populated.n1.interval <- compactLowRepeatSeg(populated.n1.interval)
  populated.n1.interval <- annotateSeg(populated.n1.interval,
                                       loc.headers=c("chr", "start.pos", "end.pos"))
  populated.n1.interval$length <- with(populated.n1.interval, (loc.end - loc.start) / 100000)
  colnames(populated.n1.interval)[4] <- 'length_(mb)'
  het.genes <- mapGenes(populated.n1.interval)
  write.table(populated.n1.interval[,c("chrom", "loc.start", "loc.end", "PNETs_(n=17)", 
                                       names(disease.list[['otb-pnet']]), pnet.ids, 'gene')],
              file="/Users/rquevedo/Desktop/pnet_heterozygous_regions.n1collapse.txt", 
              quote = FALSE, sep="\t", col.names=TRUE, row.names=FALSE, na = "")
  write.table(populated.interval[which(populated.interval$`PNETs_(n=17)` <= 1),
                                 c("chrom", "loc.start", "loc.end", "PNETs_(n=17)", 
                                   names(disease.list[['otb-pnet']]), pnet.ids, 'gene')],
              file="/Users/rquevedo/Desktop/pnet_heterozygous_regions.n1raw.txt", 
              quote = FALSE, sep="\t", col.names=TRUE, row.names=FALSE, na = "")
  
  populated.n.high.interval <- populated.interval[which(populated.interval$`PNETs_(n=17)` >= 13),]
  populated.n.high.interval <- compactLowRepeatSeg(populated.n.high.interval)
  populated.n.high.interval <- annotateSeg(populated.n.high.interval,
                                       loc.headers=c("chr", "start.pos", "end.pos"))
  populated.n.high.interval$length <- with(populated.n.high.interval, 
                                           (loc.end - loc.start) / 100000)
  colnames(populated.n.high.interval)[4] <- 'length_(mb)'
  loh.genes <- mapGenes(populated.n.high.interval)
  write.table(populated.n.high.interval[,c("chrom", "loc.start", "loc.end", "PNETs_(n=17)", 
                                       names(disease.list[['otb-pnet']]), pnet.ids, 'gene')],
              file="/Users/rquevedo/Desktop/pnet_heterozygous_regions.n.high.collapse.txt", 
              quote = FALSE, sep="\t", col.names=TRUE, row.names=FALSE, na = "")
  
  save(het.genes, loh.genes, file="/Users/rquevedo/Desktop/NETseq_loh-het-genes.Rdata")
  
  
  ### Formats the netseq.df to append to the main cgh.df
  netseq.df <- as.data.frame(matrix(ncol=10, 
                                    nrow=length(list.files(path=net.seq.dir, pattern="T_2-N.gz_segments.txt"))), 
                             stringsAsFactors=FALSE, check.names=FALSE)
  colnames(netseq.df) <- colnames(cgh.df)
  rownames(netseq.df) <- gsub("-T.+$", "", list.files(path=net.seq.dir, pattern="T_2-N.gz_segments.txt"))
  netseq.df$Case <- NA
  netseq.df$Sex <- NA
  netseq.df$Age <- NA
  netseq.df$NET <- c("PNET", "GINET", "PNET", "GINET", "GINET", "PNET", "PNET")
  netseq.df$F.stat <- "Unk"
  netseq.df$Met <- "Liver"
  netseq.df$Losses <- loss.list[rownames(netseq.df)]
  netseq.df$Gains <- gain.list[rownames(netseq.df)]
  netseq.df$Met.stat <- "Yes"
  netseq.df$Study.id <- "NETSEQ"
  cgh.df <- rbind(cgh.df, netseq.df)
  
  ### Formats the netseq.df to append to the main cgh.df  # - EXTENDED COHORT
  netseq.df <- as.data.frame(matrix(ncol=10,
                                    nrow=length(disease.list[['otb-pnet']])),
                             stringsAsFactors=FALSE, check.names=FALSE)
  colnames(netseq.df) <- colnames(cgh.df)
  rownames(netseq.df) <- gsub("-2-0", "-1", names(disease.list[['otb-pnet']]))
  netseq.df$Case <- NA
  netseq.df$Sex <- NA
  netseq.df$Age <- NA
  netseq.df$NET <- rep("PNET", length(disease.list[['otb-pnet']]))
  netseq.df$F.stat <- "Unk"
  netseq.df$Met <- "Liver"
  netseq.df$Losses <- loss.list[names(disease.list[['otb-pnet']])]
  netseq.df$Gains <- gain.list[names(disease.list[['otb-pnet']])]
  netseq.df$Met.stat <- c(rep("Yes", 3), rep("No", 1), rep("Yes", 2), rep("No", 1),
                          rep("Yes", 1), rep("No", 5)) #, rep("Yes", 1), rep("No", 1), rep("Yes", 8)
  netseq.df$Study.id <- "NETSEQ.ext"
  cgh.df <- rbind(cgh.df, netseq.df)

  cgh.df <- cgh.df[-which(cgh.df$NET != 'PNET'),]
}



### -------------------------------- Clustering -------------------------------- ###
### Generates concordances and unsupervised clusters between CGH profile of samples
### ---------------------------------------------------------------------------- ###
# Prepares the UCSC chrom info file 
print(paste(Sys.time(), ": Reading in UCSC Chrom Info...", sep=""))
ucsc.chrom.df <- read.csv(ucsc.chromInfo, header=T, sep="\t", check.names=FALSE)
colnames(ucsc.chrom.df)[1] <- "chrom"
ucsc.chrom.df$chrom <- gsub("^chr", "", ucsc.chrom.df$chrom)
global.bin.size <- 100000 #http://dgv.tcag.ca/dgv/app/statistics?ref=

split.cgh.list <- split(cgh.df, cgh.df$Study.id)

# Get all sample level Losses and gains in a dataframe per sample (each sample is an element in the list)
split.cgh.loss <- list()
split.cgh.gain <- list()
for(study.id in names(split.cgh.list)){
  print(study.id)
  clearNull <- function(x){
    null.idx <- tryCatch({null.idx <- which(sapply(x, is.null))},
             error=function(e){
               null.idx <- c()
             })
    if(length(null.idx) > 0){
      x1 <- x[-null.idx]
      x <- x1
    }
    return(x)
  }
  
  split.cgh.loss <- append(split.cgh.loss, lapply(clearNull(split.cgh.list[[study.id]]$Losses), function(x) getSampleDf(x,-1)))
  split.cgh.gain <- append(split.cgh.gain, lapply(clearNull(split.cgh.list[[study.id]]$Gains), function(x) getSampleDf(x,1)))
}

# Combine the Gains and Losses together into a single dataframe per element/sample in a list
split.cgh <- appendList(split.cgh.loss, split.cgh.gain)
global.chrom <- c(1:22, "X", "Y")

# Order the dataframes by chromosomes and start positions
split.cgh <- lapply(split.cgh, function(x){
  x$Chromosome <- factor(x$Chromosome, global.chrom)
  x[with(x, order(Chromosome)),]
} )

#Splits the entire genome into "global.bin.size" sized bins and populates them with Gains (+1), losses (-1) or neutral (0)
chrom.bins.list <- generateGenomeBins(ucsc.chrom.df ,global.bin.size)
split.cgh.bins <- lapply(split.cgh, function(x) populateGenomeBins(x, chrom.bins.list, "stat"))
names(split.cgh.bins) <- gsub("-2-0", "-1", names(split.cgh.bins))

#Get % of Genomic-unstable 
ci.frac <- unlist(lapply(split.cgh.bins, function(x) getCIFraction(x, global.bin.size)))
ci.ord <- match(rownames(cgh.df), names(ci.frac))
cgh.df$ci.frac <- as.numeric(ci.frac[ci.ord])


# Find the concordance between all CGH profiles using a perfect match
cgh.conc <- lapply(split.cgh.bins, function(a) {
  a.stat <- a[,'stat']
  lapply(split.cgh.bins, function(b){
    b.stat <- b[,'stat']
    # Finds all matches between 1=1, -1=-1, and 0=0
    len.ab <- length(which(a.stat == b.stat))
    # Finds all the 0=0 matches to remove from counts
    zero.ab <- length(which(a.stat == b.stat & a.stat == 0))
    if(length(a.stat) == zero.ab){
      len.ab <- 0
    } else {
      #len.ab <- (len.ab - zero.ab)/(length(a.stat) - zero.ab)
      len.ab <- len.ab / length(a.stat)
    }
    return(len.ab)
  }) 
})
cgh.conc.mat <- do.call("rbind", cgh.conc)



### -------------------------------- Ordering -------------------------------- ###
###Order the CGH dataframe
### -------------------------------------------------------------------------- ###
#cgh.df.bkup <- cgh.df
#cgh.df <- cgh.df[order(cgh.df$ci.stat),]
x <- hclust(dist(cgh.conc.mat))
cgh.df <- cgh.df[match(x$labels[x$order], rownames(cgh.df)),]

# Assigns group.ids to different samples in the CGH dataframe
k.clust <- 5  # 7
hclust.groups <- cutree(x, k=k.clust)
hclust.ord <- match(rownames(cgh.df), names(hclust.groups))
cgh.df$hclust <- as.integer(hclust.groups[hclust.ord])

## Hardcoded data to visualize CI Fraction for each group
if(1==0){
  require(scales)
  clust.split <- split(cgh.df, cgh.df$hclust)
  clust.col <- rainbow(length(clust.split))
  pdf("ci_frac.clusters.pdf")
  plot(0, type="n", ylim=c(0,5), xlim=c(-0.1,1),
       ylab="log2 counts", xlab="CI-fractions", main="Genomic fraction of chromosomal instability")
  dens.list <- list()
  for(each.cgh.clust in c(1:length(clust.split))){  #length(clust.split)
    temp.dens <- density(clust.split[[each.cgh.clust]]$ci.frac, from=-0.1, to=1, n=1000)
    lines(temp.dens$x, (log2(temp.dens$y)+1), col=clust.col[each.cgh.clust])
    dens.list[[each.cgh.clust]] <- temp.dens
  }
  int.points <- lapply(dens.list, function(d1){
                                              lapply(dens.list, function(d2){ 
                                                                             poi <- which(diff(log2(d1$y) < log2(d2$y)) != 0) 
                                                                             d1$x[poi]
                                                                             })
                                              })
  int.points <- lapply(int.points, function(x) lapply(x, function(y) y[which(y > 0.05 & y < 0.5)]))
  int.points <- unlist(lapply(int.points, function(x) min(unlist(x))))
  int.points <- unique(int.points)
  abline(v = int.points, col="grey", lty='dashed')
  abline(v = mean(int.points), col="black")
  text(x=(mean(int.points) + 0.01), y=3, labels=(round(mean(int.points), 3)), adj=0)
  
  legend(0.8,4, names(clust.split),
         lty=c(1,1), lwd=c(2.5,2.5),col=clust.col)
  dev.off()
}

#Creates a contigency table and tests for significance between High/Low-CI and metadata categories
# mean(int.points) = 0.117
require(mclust)
fit <- Mclust(cgh.df$ci.frac)
if(fit$G > 2) warning("Number of EM clusters is greater than 2!")
max.idx <- which(fit$parameters$mean >= 0.2)
high.idx <- which(fit$classification %in% names(max.idx))
cgh.df$ci.stat <- "Low-CI"
cgh.df$ci.stat[high.idx] <- 'High-CI'
#cgh.df$ci.stat <- ifelse(cgh.df$ci.frac >= 0.117, "High-CI", "Low-CI")
hclusters <- cgh.df$ci.stat
metadata.clust <- c("Study.id", "NET", "F.stat", "Met.stat")
assoc.df <- data.frame(matrix(ncol=6, nrow=0))
for(each.clust in metadata.clust){
  y <- assocstats(table(cgh.df[,each.clust], hclusters,
                   dnn=c(each.clust, "GI.stat")))
  assoc.st <- data.frame(t(c("CI fraction", each.clust, 
                             round(y$chisq_tests['Pearson',][1], 2), round(y$chisq_tests['Pearson',][2], 1), round(y$chisq_tests['Pearson',][3], 8),
                             "cramer"=round(y$cramer,4))))
  assoc.df <- rbind(assoc.df, assoc.st)
}
# Hardcoded colnames and relabelling of category-2 for tables
colnames(assoc.df) <- c("Category-1", "Category-2", "X^2", "d.f.", "Pr(>X^2)", "V")
assoc.df[,'Category-2'] <- c("Study ID", "NET type", "Functional status", "Metastasis status")

if(0==1){
  # Hardcoded generation of latex tables for the assoc.df and significant associations [ study.id and met.stat ]
  #X\textsuperscript{2}
  print(xtable(assoc.df, include.rownames=FALSE ))
  
  print(xtable(table(cgh.df[,'Study.id'], hclusters,
                     dnn=c('Study ID', "CI Fraction"))))
  
  print(xtable(table(cgh.df[,'Met.stat'], hclusters,
                     dnn=c('Study ID', "CI Fraction"))))
  
  print(xtable(table(cgh.df[,'F.stat'], hclusters,
                     dnn=c('Functionality', "CI Fraction"))))
  
  print(xtable(table(cgh.df[,'NET'], hclusters,
                     dnn=c('NET Type', "CI Fraction"))))
}


### -------------------------------- Plotting -------------------------------- ###
### Splits the screen dynamically according to samples and plots the CGH data
### -------------------------------------------------------------------------- ###

### Plotting pre-processing and setting up screen.split matrices for samples
#Groups to split cgh.df by:
#    "Study.id"    "NET"      "F.stat"   "Met"  "ci.stat"
groups <- c("Study"="Study.id",
            "NET type"="NET", 
            "Functionality"="F.stat", 
            "Metases"="Met.stat",
            "CI-Fraction"="ci.stat")
num.of.groups <- length(groups)

# Set-up Chr-plotting split.screen
lspacer <- 0.1
rspacer <- 0.9
dim.sep <- seq(lspacer, rspacer, by=((rspacer-lspacer)/dim(cgh.df)[1]))

# Set up samples split matrix
all.samples.mat <- matrix(c(0, dim.sep,
                            lspacer, (dim.sep[-1]), 1.0,
                            0,rep(0, dim(cgh.df)[1]),0,
                            1,rep(1,dim(cgh.df)[1]), 1), 
                          ncol=4)
hclust.split <- matrix(c(0,lspacer, 0, 1,
                         lspacer, rspacer, 0, 1,
                         rspacer, 1, 0, 1), ncol=4, byrow=TRUE)

pdf("cghPlots.cistat.netseq.pdf")
# Plots the heading track above chr-plots
split.screen(matrix(c(rep(0, num.of.groups+4),
                      rep(1, num.of.groups+4),
                      c(0.000, 0.050, 0.075, 0.100, 0.125, 0.150, 0.175, 0.2, 0.9),
                      c(0.050, 0.075, 0.100, 0.125, 0.150, 0.175, 0.200, 0.9, 1.0)),
                    ncol=4))
### Hclust Plots:
hclust.screen <- num.of.groups + 4
screen(hclust.screen)
screen.id <- split.screen(hclust.split)
screen(screen.id[2])
par(mar=c(0,0,0.3,0))
cgh.hc <- hclust(dist(cgh.conc.mat))
plot(cgh.hc, labels=FALSE, hang=-1, 
     axes=FALSE, ylab="", xlab="", main="", sub="", yaxs="i", xaxs="i")
cgh.col <- brewer.pal(k.clust,"Set1")
for(i in c(1:length(cgh.df$hclust))){
  rect(xleft=(i-1), ybottom = -1, xright = i, ytop = 0, 
       col=cgh.col[cgh.df$hclust[i]], lty=0)
}
#rect.hclust(cgh.hc, k=k.clust, border="red")

### Chromosome CNA plot per sample
chr.screen <- num.of.groups+3   # num.of.groups + upper spacer + chr.space-base
screen(chr.screen)
screen.id <- split.screen(all.samples.mat)
base.screen <- min(screen.id)   # Which screen to start plotting on
for(each.sample in c(1:dim(cgh.df)[1])){
  screen(each.sample + base.screen) 
  par(mar=c(0,0,0,0))
  plot(0, type='n', 
       xlim=c(0,1), ylim=c(-t.genome.size, 0), 
       axes=FALSE,
       ylab='', xlab='', yaxs='i', xaxs='i')
  plotCNA(cna.info=cgh.df[each.sample,'Losses'], cna.col="blue", ucsc.chrom=ucsc.chrom, t.genome.size=t.genome.size)
  plotCNA(cna.info=cgh.df[each.sample,'Gains'], cna.col="red", ucsc.chrom=ucsc.chrom, t.genome.size=t.genome.size)
  abline(h = -(ucsc.chrom$chr.st.size), col="grey")
}
#Label the chromosomes
screen(base.screen)
par(mar=c(0,3.5,0,0))
plot(0, type='n', 
     xlim=c(0,1), ylim=c(-t.genome.size, 0), 
     axes=FALSE, ylab='', xlab='', yaxs='i', xaxs='i')
axis(side = 2, at = c(-t.genome.size, -(rev(ucsc.chrom$chr.st.size))), labels = rep("",23), las=2)
axis(side = 2, at = -(rev(ucsc.chrom$mid.st)), labels = c(22:1), las=2, lty=0, cex.axis=0.7)



### Screens for the colour of group-level information
#group.col <- c("spiral.morning", "qual", "survival", "old.qual2", "old.qual2")
for(each.group in c(1:num.of.groups)){
  group.screen <- each.group + 1  # +1 because of upper spacer
  screen(group.screen)
  screen.id <- split.screen(all.samples.mat)
  base.screen <- min(screen.id)    # Which screen to start plotting on
  
  
  sample.col <- convertFactorToColor(cgh.df[, groups[each.group]], groups[each.group])
  for(each.sample in c(1:length(sample.col))){
    screen(each.sample + base.screen) 
    par(mar=c(0,0,0,0))
    plot(0, type="n", xlim=c(0,1), ylim=c(0,1), axes=FALSE)
    rect(0,0,1,1, col=sample.col[each.sample], border='white')
  }
  
  #Label the groups
  screen(base.screen)
  par(mar=c(0,3.5,0,0))
  plot(0, type='n', 
       xlim=c(0,1), ylim=c(0,1), 
       axes=FALSE, ylab='', xlab='')
  axis(side = 2, at = 0.5, labels = names(groups)[each.group], las=2, lty=0, cex.axis=0.5, tck=0)
}
close.screen(all.screens=TRUE)

### Plot the legends    [hardcoded]
split.screen(c(2,3))
for(each.group in c(1:num.of.groups)){
  screen(each.group)
  sample.col.list <- convertFactorToColor(cgh.df[, groups[each.group]], groups[each.group], ret.ord=1)
    
  ord.col <-  as.vector(sample.col.list$color)
  ord.id <- as.vector(sample.col.list$id)
  if(names(groups)[each.group] == 'Functionality'){
    ord.id[which(ord.id == "+")] <- "Functional"
    ord.id[which(ord.id == "-")]  <- "Non-functional"
    ord.id[which(ord.id == "Unk")] <- "Unknown"
  }
  
  par(mar=c(0,0,0,0))
  plot(0, type="n", xlim=c(0,100), ylim=c(0,100), axes=FALSE, ylab="", xlab="")
  text(x=10, y=95, labels = '', adj = 0)
  legend(10,90,  
         ord.id,
         cex=0.7,
         fill=ord.col,
         title=names(groups)[each.group],
         bty="n", xjust=0)
}
#Add in clusters legend:
screen(num.of.groups + 1)
ord.col <- cgh.col
ord.id <- c(1:k.clust)
par(mar=c(0,0,0,0))
plot(0, type="n", xlim=c(0,100), ylim=c(0,100), axes=FALSE, ylab="", xlab="")
text(x=10, y=95, labels = '', adj = 0)
legend(10,90,  
       ord.id,
       cex=0.7,
       fill=ord.col,
       title="Clusters",
       bty="n", xjust=0)

close.screen(all.screens=TRUE)

dev.off()



# setwd('~/git/loh_cn_visualize/cgh_plotter')
# save(ucsc.chrom, hg19.cytoband, cgh.list, cgh.df, split.cgh.bins, file="cghRawData.netseq.Rdata")


###### Plot the densities of the different clusters and find the minimum intersection:
pdf("lohCIdensity.pdf")
split.screen(c(2,1))
screen(1)
cgh.clust.list <- split(cgh.df, cgh.df$hclust)
#cgh.clust.list <- cgh.clust.list[-5]
cgh.clust.dens <- lapply(cgh.clust.list, function(x) density(x$ci.frac, from=-0.1, to=1))
plot(0,0, ylim=c(0,1.5), xlim=c(0,1), type='n', yaxt='n',
     xlab="LOH Genomic Fraction", ylab="Density")
for(i in names(cgh.clust.dens)){
  lines(x=cgh.clust.dens[[i]]$x, y=log10(cgh.clust.dens[[i]]$y + 1), col=cgh.col[as.integer(i)])
}

all.vals <- c()
for(a in names(cgh.clust.dens)){
  x <- cgh.clust.dens[[a]]
  tpoi <- c(200)
  for(b in names(cgh.clust.dens)){
    y <- cgh.clust.dens[[b]]
    poi <- which(diff(x$y < y$y) != 0)
    poi <- poi[which((x$x[poi] < 0.3) & (x$x[poi] > 0.07))]
    tpoi <- c(tpoi, poi)
  }
  if(a %in% '5'){
    tpoi <- tpoi[-match(min(tpoi), tpoi)]
  }
  text(x=(x$x[min(tpoi)] + 0.03), 
       y=(log10(x$y[min(tpoi)] + 1) + 0.05),
       labels=round(x$x[min(tpoi)], 3), cex=0.3)
  points(x=x$x[min(tpoi)], y=log10(x$y[min(tpoi)] + 1))
  all.vals <- c(all.vals, x$x[min(tpoi)])
}
abline(v=mean(all.vals))
text(x=(mean(all.vals) + 0.075), y=1.2, labels=paste("Lower GI: ", round(mean(all.vals), 3), sep=""), cex=0.6)

rect(xleft = mean(all.vals), ybottom = 1.4, xright = 1, ytop = 1.5, col=alpha("red", 0.5))
rect(xleft = 0, ybottom = 1.4, xright = mean(all.vals), ytop = 1.5, col=alpha("blue", 0.5))
text(x= 0.6, y=1.45, labels="High-CI", cex=0.6)
text(x= 0.05, y=1.45, labels="Low-CI", cex=0.6)

axis(2, at=seq(0, 1.5, by=0.5),
     labels=c(0, (round(10^(seq(0.5, 1.5, by=0.5)), 1) + 1)),   las=2)
close.screen(all.screens=TRUE)
dev.off()

save(cgh.df, file="~/git/loh_cn_visualize/data/cgh_df.RData")
