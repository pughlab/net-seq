library(scales)
library(mixdist)
library(intervals)
library(VariantAnnotation)
library(rowr)
library(metafor)

source('~/git/net-seq/af_plotter/src/visualizationAf.R')
load('~/git/net-seq/af_plotter/src/data/afPlotter.ref.Rdata') # Contains: clinvar.df, refgene.df, chrom.df

###############
##### Set-up
global.bin.size <- 500000  # genomic bin size
global.chrom <- c(1:22) # chromosome ID
assembly <- 'GRCh37'  # GRCh37 or GRCh38
source('~/git/net-seq/af_plotter/src/snpProcess.R')
source('~/git/net-seq/af_plotter/src/binProcess.R')
source('~/git/net-seq/af_plotter/src/misc.R')
source('~/git/net-seq/shallow_wgs/src/physicalCov.R')
source('~/git/net-seq/shallow_wgs/src/misc.R')
source('~/git/net-seq/shallow_wgs/src/ORanalysis.R')
load('~/git/net-seq/shallow_wgs/data/referenceDist.Rdata')

ecdf.metric <- "median" # mean or median
ecdf.thresh <- 0.25
ecdf.stroma.thresh <- TRUE  # If true, will use a threshold corresponding to stromal content
collapse.type <- 'seg'  # seg or bin

# Optional: Odds Ratio Test
doOddsRatioTest <- TRUE
generateSequenzaFiles <- TRUE
plotStromaLohRegression <- TRUE
if(doOddsRatioTest){
  anno.file <- '~/git/net-seq/shallow_wgs/data/net.wgsAnno.csv'
  anno.df <- read.csv(anno.file, header=TRUE, stringsAsFactors=FALSE)
  group <- 'met.stat'   # met.stat or net.type
  group.id <- c("met", "nonMet")  # ("met, "nonMet") or ("PNET", "GINET")
}


# Generate preliminary reference dataframes and genomic bins
clinvar.df <- formatClinvarDf(clinvar.df, assembly)
chrom.bins.list <- generateGenomeBins(ucsc.chrom.df ,global.bin.size)

setwd('~/git/net-seq/shallow_wgs/data')
all.net.id <- c("004", "001a", "001b", "002", "003","005", 
                "006", "007", "008", "009", "011", "012", "013", "014", 
                "015", "016", "017", "021", "024", "025", "026", "027", 
                "028", "029", "030", "031", "032", "033", "034", "036", "037")
pnet.id <- c('001a', '001b', '002', '003', '005', '011', '012', '013', '014', '015',
             '017', '025', '028', '029', '030', '031', '032', '034', '035')
ginet.id <- c('004', '006', '007', '008', '009', '010', '016', '018', '019', '020',
              '021', '022', '023', '024', '026', '027', '033', '036', '037')
tumour.stroma <- c('001a'=0.05, '001b'=0.10, '002'=0.20, '003'=0.50, '004'=0.40, '005'=0.75,
                   '006'=0.10, '007'=0.25, '008'=0.40, '009'=0.40, '011'=0.20, '012'=0.40,
                   '013'=0.20, '014'=0.30, '015'=0.30, '016'=0.20, '017'=0.30, '021'=0.20,
                   '024'=0.30, '025'=0.40, '026'=0.50, '027'=0.20, '028'=0.40, '029'=0.30,
                   '030'=0.25, '031'=0.40, '032'=0.30, '033'=0.50, '034'=0.50, '036'=0.20, '037'=0.20)

#Samples to remove upon review from our resident neuroendocrine pathologist Dr. Sylvia Asa
pathology.blacklist <- paste("0", c("01a", "01b", "08", "12", "14", "26", "28", "36"), sep="") 
all.net.id <- all.net.id[which(!all.net.id %in% pathology.blacklist)]

disease.list <- list('otb-pnet'=list(),
                     'otb-ginet'=list())
if(collapse.type %in% 'bin') print("Disease DF based on 500kB bin-level resolution") else 
  print("Disease DF based on collapsed CN segment-level resolution")
sample.values <- list()

for(net.id in all.net.id){
  if(ecdf.stroma.thresh){
    ecdf.thresh <- as.numeric(tumour.stroma[net.id])
  }
  file.list <- c(wgs_t=paste('NET-2-', net.id, '.snpOut.vcf', sep=""),
                 wgs_n='allNET.snpOut.vcf',
                 wgs_tcov=paste('NET-2-', net.id, '.bedpe.cov', sep=""))
  
  
  
  
  ###############
  ##### Main
  
  ### 1: Process the physical Coverage:
  wgsCov.df <- read.table(file.list['wgs_tcov'], header=FALSE, sep="\t", stringsAsFactors=FALSE)
  
  colnames(wgsCov.df) <- c("chrom", "start.pos", "end.pos", "count")
  wgsCov.df <- wgsCov.df[with(wgsCov.df, order(chrom, start.pos)),]
  wgsCov.list <- split(wgsCov.df, wgsCov.df$chrom)
  
  # Obtain the mixed Gaussian EM estimated mu's
  fitpro <- getEMestimates(wgsCov.df, 200)
  
  # Fit observed tumour counts to all binomial models centered around EM mu's
  mu.model <- fitpro$parameters$mu
  mu.model[which(mu.model < 0)] <- 0
  prob.fit <- lapply(mu.model, function(mu.cnt){
    apply(wgsCov.df, 1, dbinomCounts, success.cnt=mu.cnt, bsize=global.bin.size)
  })
  prob.fit.df <- do.call("cbind", prob.fit)
  
  # Assign the proper copy-number to each peak:
  #   Default is a sequential order  (e.g. colnames(prob.fit.df) <- c("a", "b", "c", "d"))
  fitted.cn <- fitCnLabels(fitpro$parameters)
  colnames(prob.fit.df) <- fitted.cn$abs.cn
  #colnames(prob.fit.df) <- as.character(c(0:(length(mu.model)-1)))
  
  # Extract the total-CN based on the highest probability for the fit mu-binomial models
  prob.fit.cn <- apply(prob.fit.df, 1, getMaxProb, colnames(prob.fit.df))
  prob.fit.cn.df <- do.call("rbind", prob.fit.cn)[,1, drop=FALSE]
  colnames(prob.fit.cn.df) <- c("copy.number")
  wgsCov.df <- cbind(wgsCov.df, prob.fit.df, prob.fit.cn.df)
  
  
  
  ### 2: Process the allelic fractions: 
  # Extract tumor AF and bin
  each.vcf <- 'wgs_t'
  tumor.af.list <- generateChrList(file.list[each.vcf], getwd(), snp.caller='mutect', clinvar.stat=FALSE, col.id=2)
  tumor.af.list.bin <- binSnps(tumor.af.list, rm.hom=FALSE, shallow.wgs=TRUE)
  
  # Extract normal AF and bin
  each.vcf <- 'wgs_n'
  # dna.af.list <- generateChrList(file.list[each.vcf], getwd(), snp.caller='mutect', clinvar.stat=FALSE, col.id=2)
  # dna.af.list.bin <- binSnps(dna.af.list, rm.hom=FALSE, shallow.wgs=TRUE)
  dna.af.list.bin <- split(normal.df, normal.df$chr)
  
  # Obtain the quantile distribution of tumor to normal ECDF
  het.cdf.list <- lapply(dna.af.list.bin, function(x) ecdf(x[,'het']))
  quant.het <- lapply(tumor.af.list.bin, function(x) het.cdf.list[['chr1']](x[,'het']))
  for(each.chr in names(tumor.af.list.bin)){
    tumor.af.list.bin[[each.chr]] <- cbind(tumor.af.list.bin[[each.chr]], 
                                           quantile=quant.het[[each.chr]])
    colnames(tumor.af.list.bin[[each.chr]])[1] <- 'chrom'
  }
  
  
  ### 3: Intersect the AF and DoC into one dataframe:
  wgsCov.list <- split(wgsCov.df, wgsCov.df$chrom)
  wgsCovAf.list <- list()
  for(each.chrom in names(tumor.af.list.bin)){
    wgsCovAf.list[[each.chrom]] <- merge(x=wgsCov.list[[each.chrom]], 
                                         y=tumor.af.list.bin[[each.chrom]], 
                                         by=c("chrom", "start.pos", "end.pos"), all=TRUE)
  }
  wgsCovAf.df <- do.call("rbind", wgsCovAf.list)
  
  wgsCovAf.df$copy.number <- as.integer(as.character(wgsCovAf.df$copy.number))
  wgsCovAf.seg <- collapseCov(wgsCovAf.df)
  
  wgsCovAfDf.list <- wgsCovAf.list
  wgsCovAfSeg.list <- split(wgsCovAf.seg, wgsCovAf.seg$chrom)
  
  chr.values <- lapply(wgsCovAfSeg.list, function(x){
    med.val <- median(rep(x$af.quant.med, x$num.of.bins), na.rm=TRUE)
    mean.val <- mean(rep(x$af.quant.mu, x$num.of.bins), na.rm=TRUE)
    return(c("med"=med.val, "mu"=mean.val))
  })
  for(each.chr in names(chr.values)){
    sample.values[[each.chr]] <- rbind(sample.values[[each.chr]], 
                                       c(chr.values[[each.chr]], tumour.stroma[net.id]))
    rownames(sample.values[[each.chr]])[dim(sample.values[[each.chr]])[1]] <- net.id
  }
  
  
  ### 4: Assemble contigency tables for each sample and each chromosome
  if(doOddsRatioTest){
    # Get the representative AF-percentile state across the entire chromosome
    # Report 1 if it's below Q1 and 0 if above Q1
    if(ecdf.metric %in% 'mean') bin.m <- 'af.quant.mu' else bin.m <- 'af.quant.med'
    summ.af.percentile <- lapply(wgsCovAfSeg.list, getChromSnpNeg, bin.metric=bin.m)
    neg.df <- do.call("rbind", lapply(summ.af.percentile, function(x) if(x < ecdf.thresh) 1 else 0))
    neg.df <- neg.df[order(as.integer(gsub("chr", "", rownames(neg.df)))),, drop=FALSE]
    colnames(neg.df) <- net.id
    
    # Get the representative CN state across the entire chromosome
    summ.cn.val <- lapply(wgsCovAfSeg.list, getChromCnNeg)
    cn.df <- round(do.call("rbind", summ.cn.val),0)
    cn.df <- cn.df[order(as.integer(gsub("chr", "", rownames(cn.df)))),, drop=FALSE]
    colnames(cn.df) <- net.id
    
    if(net.id %in% all.net.id[1]){
      cn.state.df <- cn.df
      snp.deficient.df <- neg.df
    } else {
      cn.state.df <- cbind.fill(cn.state.df, cn.df, fill=NA)
      snp.deficient.df <- cbind.fill(snp.deficient.df, neg.df, fill=NA)
    }
  }
  
  
  ### 5: Create DoC and AF Plots:
  dir.create(file.path("~/Desktop/netseq/net_extended/plots"), recursive = TRUE, showWarnings = FALSE)
  png(paste("~/Desktop/netseq/net_extended/plots/shallowWgsCnv.NET", net.id, ".png", sep=""), 
      width=20, height=8, units='in', res = 300)
  split.screen(c(5,1))
  screen(2) # Coverage Screen
  doc.screen.id <- split.screen(c(1, (22 + 2)))
  doc.screen.id <- doc.screen.id[-c(1, length(doc.screen.id))]
  chr.cnt <- 1
  for(i in doc.screen.id){
    screen(i)
    par(mar=c(1,0.2,0.5, 0))
    plot(wgsCov.list[[paste("chr", chr.cnt, sep="")]]$count, 
         ylim=c(1,3000), ylab='', xaxt='n', xlab=chr.cnt,
         yaxt=if((i) == doc.screen.id[1]) 't' else 'n')
    chr.cnt <- chr.cnt + 1
  }
  
  screen(3) # AF Screen
  af.screen.id <- split.screen(c(1, (22 + 2)))
  af.screen.id <- af.screen.id[-c(1, length(af.screen.id))]
  chr.cnt <- 1
  for(i in af.screen.id){
    screen(i)
    par(mar=c(1,0.2,0.5, 0))
    if(is.null(wgsCovAfDf.list[[paste("chr", chr.cnt, sep="")]])){
      plot(0, type='n')
    } else {
      plot(wgsCovAfDf.list[[paste("chr", chr.cnt, sep="")]]$quantile, ylim=c(0,1), 
           ylab='', xaxt='n', xlab=chr.cnt,
           yaxt=if((i) == af.screen.id[1]) 't' else 'n')
    }
    
    chr.cnt <- chr.cnt + 1
  }
  
  
  
  colfunc <- colorRampPalette(c("blue", "black", "red"))(100)
  
  screen(4) # AF Screen
  seg.screen.id <- split.screen(c(1, 22 + 2))
  seg.screen.id <- seg.screen.id[-c(1, length(seg.screen.id))]
  chr.cnt <- 1
  for(i in seg.screen.id){
    screen(i)
    par(mar=c(1,0.2,0.5, 0))
    if(is.null(wgsCovAfDf.list[[paste("chr", chr.cnt, sep="")]])){
      plot(0, type="n")
    } else {
      plot(x=wgsCovAfDf.list[[paste("chr", chr.cnt, sep="")]]$end.pos,
           y=wgsCovAfDf.list[[paste("chr", chr.cnt, sep="")]]$copy.number, type='n', 
           ylim=c(0,6), ylab='', xaxt='n', xlab=chr.cnt, las=2, cex.axis=0.75,
           yaxt=if((i) == seg.screen.id[1]) 't' else 'n')
      apply(wgsCovAfSeg.list[[paste("chr", chr.cnt, sep="")]], 1, function(x){
        rect(xleft=as.integer(as.character(x['start.pos'])),
             ybottom=(as.integer(as.character(x['copy.number'])) - 0.2),
             xright=as.integer(as.character(x['end.pos'])),
             ytop=(as.integer(as.character(x['copy.number'])) + 0.2),
             col=alpha(colfunc[round(as.numeric(as.character(x['af.quant.mu'])) * 100, 0)], 0.5))
      })
    }
    chr.cnt <- chr.cnt + 1
  }
  
  screen(5) # Chr labelling screen
  seg.screen.id <- split.screen(c(1, 22 + 2))
  seg.screen.id <- seg.screen.id[-c(1, length(seg.screen.id))]
  chr.cnt <- 1
  for(i in seg.screen.id){
    screen(i)
    par(mar=c(1,0.2,0.5, 0))
    plot(0, xlim=c(0,10), ylim=c(0,10), axes=FALSE,
         type='n', xaxt='n', yaxt='n', ylab='', xlab='')
    text(x = 5, y=9.5, chr.cnt)
    chr.cnt <- chr.cnt + 1
  }
  
  close.screen(all.screens=TRUE)
  dev.off()
  png(paste("~/Desktop/netseq/net_extended/plots/shallowWgsCnv.NET", net.id, ".fitpro.png", sep=""), 
      width=20, height=8, units='in', res = 300)
  plot(fitpro)
  dev.off()
  
  write.table(fitted.cn, file=paste("~/Desktop/netseq/net_extended/plots/shallowWgsCnv.NET", net.id, ".txt", sep=""),
              sep="\t", quote=FALSE, col.names=TRUE, row.names=FALSE)
  ### 5: Assemble disease.list for plotting on aggregateLoh.R
  dl.net.id <- paste("NET-2-", net.id, sep="")
  if(net.id %in% pnet.id) net.type <- 'otb-pnet' else net.type <- 'otb-ginet'
  print(paste(dl.net.id, "; ECDF Threshold = ", ecdf.thresh, sep=""))
  if(collapse.type %in% 'seg'){
    diseaselist.df <- do.call("rbind", lapply(wgsCovAfSeg.list, genDiseaseDF, type='seg'))
    disease.list[[net.type]][[dl.net.id]]  <- as.data.frame(diseaselist.df)
  } else {
    diseaselist.df <- do.call("rbind", lapply(wgsCovAfDf.list, genDiseaseDF, type='bin'))
    disease.list[[net.type]][[dl.net.id]] <-  collapseCopyState(diseaselist.df)
  }
}

save(disease.list, file=paste('~/Desktop/netseq/net_extended/plots/plots/otbData.aggregateLoh.', 
                         ecdf.metric, 'LOH.Rdata', sep=""))


if(doOddsRatioTest){
  colnames(snp.deficient.df) <- paste("X", all.net.id, sep="")
  
  ORgroups <- list("pnetMet"=list("group" = c('net.type', 'met.stat'),
                                  "group.id" = c("PNET.met", "PNET.nonMet")),
                   "metStat"=list(group = 'met.stat',  
                                  group.id = c("met", "nonMet")),
                   "netStat"=list(group = 'net.type',
                                  group.id = c("PNET", "GINET")))
  
  lapply(ORgroups, genORTable)
}

if(generateSequenzaFiles){
  outdir <- '~/Desktop/net_OTB'
  
  for(each.net in c("pnet", "ginet")){
    dir.create(file.path(outdir, paste(each.net, ecdf.metric, sep="-")), showWarnings = FALSE)
    
    for(each.sample in names(disease.list[[paste("otb", each.net, sep="-")]])){
      sample.df <- disease.list[[paste("otb", each.net, sep="-")]][[each.sample]]
      sample.df$length <- as.integer(as.character(sample.df$End.bp)) - as.integer(as.character(sample.df$Start.bp))
      sample.df$Sample <- rep(each.sample, nrow(sample.df))
      sample.df <- sample.df[,c("Sample", "Chromosome", "Start.bp", "End.bp", "modal_A1", "modal_A2", "LOH", "length")]
      sample.df <- sample.df[order(as.integer(gsub("chr", "", sample.df$Chromosome))),]
      write.table(sample.df, file=file.path(outdir, paste(each.net, ecdf.metric, sep="-"), 
                                            paste("shallowWgsCnv.", gsub("-2-", "", each.sample), ".sequenza.txt", sep="")),
                  quote=FALSE, col.names=TRUE, row.names=FALSE, sep="\t")
    }
  }
}

if(plotStromaLohRegression){
  pdf("~/Desktop/netseq/net_extended/plots/regressionPlots.pdf")
  split.screen(c(11,2))  
  screen.cnt <- 0
  for(each.chrom in paste("chr", c(1:22), sep="")){
    screen.cnt <- screen.cnt + 1
    screen(screen.cnt)
    par(mar=c(0.2,0.25,1,0.25))
    
    stroma.purity.df <- as.data.frame(sample.values[[each.chrom]])
    colnames(stroma.purity.df) <- c("med", "mu", "stroma")
    pnet.rows <- which(rownames(stroma.purity.df) %in% pnet.id)
    stroma.purity.df.tmp <- stroma.purity.df[pnet.rows,]
    plot(med ~ stroma, data=stroma.purity.df.tmp, xlim=c(0,1), ylim=c(0,1), 
         xaxt='n', yaxt='n', main=each.chrom, cex=0.5, col="red")
    abline(lm(med ~stroma, data=stroma.purity.df.tmp), col="red")
    #text(x=stroma.purity.df.tmp$stroma, y=stroma.purity.df.tmp$med, labels=rownames(stroma.purity.df.tmp), col="red")
    
    ginet.rows <- which(rownames(stroma.purity.df) %in% ginet.id)
    stroma.purity.df.tmp <- stroma.purity.df[ginet.rows,]
    points(med ~ stroma, data=stroma.purity.df.tmp, 
           xlim=c(0,1), ylim=c(0,1), cex=0.5, col="blue")
    abline(lm(med ~stroma, data=stroma.purity.df.tmp), col="blue")
    #text(x=stroma.purity.df.tmp$stroma, y=stroma.purity.df.tmp$med, labels=rownames(stroma.purity.df.tmp), col="blue")
    
  }
  close.screen(all.screens=TRUE)
  dev.off()
}
