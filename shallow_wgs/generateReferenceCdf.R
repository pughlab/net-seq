library(scales)
library(mixdist)
library(intervals)
library(VariantAnnotation)
library(rowr)

source('~/git/net-seq/af_plotter/visualizationAf.R')
load('~/git/net-seq/af_plotter/data/afPlotter.ref.Rdata') # Contains: clinvar.df, refgene.df, chrom.df

###############
##### Set-up
global.bin.size <- 500000  # genomic bin size
global.chrom <- c(1:22) # chromosome ID
assembly <- 'GRCh37'  # GRCh37 or GRCh38
source('~/git/net-seq/af_plotter/snpProcess.R')
source('~/git/net-seq/af_plotter/binProcess.R')
source('~/git/net-seq/af_plotter/misc.R')
source('~/git/net-seq/shallow_wgs/src/physicalCov.R')
source('~/git/net-seq/shallow_wgs/src/misc.R')

# Generate preliminary reference dataframes and genomic bins
clinvar.df <- formatClinvarDf(clinvar.df, assembly)
chrom.bins.list <- generateGenomeBins(ucsc.chrom.df ,global.bin.size)


setwd('/Users/rquevedo/Desktop/PughLab/NET-seq/otb_samples/net_OTB/data')
normal.df <- data.frame()
all.net.id <- c("001a","001b","002","003","004","005","006","007","008","009",
                "011","012","013","014","015","016","017","021","024","025",
                "026","027","028","029","030","031","032","033","034","036","037")

#Samples to remove upon review from our resident neuroendocrine pathologist Dr. Sylvia Asa
pathology.blacklist <- paste("0", c("01a", "01b", "08", "12", "14", "26", "28", "36"), sep="") 
all.net.id <- all.net.id[which(!all.net.id %in% pathology.blacklist)]

for(net.id in all.net.id){
  file.list <- c(wgs_n=paste('NET-2-', net.id, '.snpOut.vcf', sep=""))
  
  normal.af.list <- generateChrList(file.list['wgs_n'], getwd(), snp.caller='mutect', clinvar.stat=FALSE, col.id=2)
  normal.af.list.bin <- binSnps(normal.af.list, rm.hom=FALSE, shallow.wgs=TRUE)
  normal.af.df <- do.call("rbind", normal.af.list.bin)
  
  if(net.id %in% '004'){
    normal.df <- normal.af.df[,c("chr", "start.pos", "end.pos", "het")]
  } else {
    normal.df <- cbind.fill(normal.df, normal.af.df[,c("het"), drop=FALSE], fill=NA)
  }
  colnames(normal.df)[dim(normal.df)[2]] <- net.id
}

colnames(normal.df) <- gsub("^X", "", colnames(normal.df))
avg.het.cnt <- apply(normal.df, 1, function(x) median(as.numeric(as.character(x[all.net.id])), na.rm=TRUE))
normal.df <- cbind(normal.df, "het"=avg.het.cnt)
save(normal.df, file="~/git/net-seq/shallow_wgs/data/referenceDist.Rdata")


