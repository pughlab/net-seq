####################################
######  formatSqzaToAbs.R
######    Purpose: Reformats the Sequenza style _segments.txt file
######      to give the Absolute style segment files for the purpose
######      of plotting LOH in the pancancer_loh.R script and other
######      related LOH-suite R scripts
######    Arguments:
######      [1] - segment file (i.e. net-001-t_2_segments.txt)
######      [2] - sample name (i.e. NET-001-T_2)
######      [3] - type (default=sequenza, afplotter)
######      [4] - cell line ploidy (1,2,3,4)
####################################


vis <- 0
#net.f <- '/Users/rquevedo/Desktop/purity_tables/dna/purity_70_BON1_baseline_ploidy3.vcf.gz.txt'
type <- 'sequenza'
args <- commandArgs(trailingOnly = TRUE)
net.f <- args[1]
net.df <- read.csv(net.f, sep="\t", quote="\"")
net.df <- cbind(Sample=rep(args[2], dim(net.df)[1]), net.df)
if(length(args) > 2){
  type <- args[3]
}


if(type == 'sequenza'){
  #Reheader sequenza to absolute standards
  colnames(net.df) <- gsub("chromosome", "Chromosome", colnames(net.df))
  colnames(net.df) <- gsub("start.pos", "Start.bp", colnames(net.df))
  colnames(net.df) <- gsub("end.pos", "End.bp", colnames(net.df))
  colnames(net.df) <- gsub("A", "modal_A1", colnames(net.df))
  colnames(net.df) <- gsub("B", "modal_A2", colnames(net.df))
  
} else if (type == 'afplotter'){
  colnames(net.df) <- gsub("chrom", "Chromosome", colnames(net.df))
  colnames(net.df) <- gsub("loc.start", "Start.bp", colnames(net.df))
  colnames(net.df) <- gsub("loc.end", "End.bp", colnames(net.df))
  
  ploidy <- args[4]
  net.df <- cbind(net.df, modal_A1=round(net.df$seg.mean.a), 
                          modal_A2=as.integer(gsub(1, ploidy, round(net.df$seg.mean.b))))
}
# Add the LOH column
loh <- apply(net.df, 1, function(x) if(x['modal_A1'] == 0 || x['modal_A2'] == 0) 1 else  0 )
#loh <- rep(0, dim(net.df)[1])
net.df <- cbind(net.df, LOH=loh)

# Add the length column
length <- apply(net.df, 1, function(x) as.numeric(x['End.bp']) - as.numeric(x['Start.bp']))
net.df <- cbind(net.df, length=length)

# select for only the needed columns:
net.df <- net.df[,c("Sample", "Chromosome", "Start.bp", "End.bp", "modal_A2f", "N.modal_A2modal_A1F",
                    "sd.modal_A2modal_A1F", "depth.ratio", "N.ratio", "sd.ratio", 
                    "CNt", "modal_A1", "modal_A2", "LPP", "LOH", "length")]

write.table(net.df, file=paste(args[2], "_segments.txt", sep=""), quote=FALSE, sep = "\t", row.names = FALSE, col.names = TRUE)





############################################################
# Visualization of the LOH Segments... Rough.  Just used to verify.  Refer to other scripts for better visualizations
if(vis == 1){
  ucsc.chrom.info.file <- '/Users/rquevedo/Desktop/bhk_lab/reference/ucsc.hg19.chromInfo.txt'     #hg19
  
  ucsc.chrom <- read.csv(ucsc.chrom.info.file, header=TRUE, sep="\t", check.names=FALSE)
  colnames(ucsc.chrom) <- gsub("#", "", colnames(ucsc.chrom))
  ucsc.chrom <-  ucsc.chrom[regexpr("_", ucsc.chrom$chrom) == -1,]
  ucsc.chrom <- ucsc.chrom[regexpr("^chrm$", ucsc.chrom$chrom, ignore.case = TRUE, perl=TRUE) == -1,]
  
  num.chr <- length(unique(net.df$Chromosome))
  net.list <- split(net.df, net.df$Chromosome)
  
  pdf(paste(args[2], "_loh.pdf", sep=""))
  for(each.chr in unique(net.df$Chromosome)){
    size <- ucsc.chrom[match(each.chr, ucsc.chrom$chrom),]$size
    plot(0, type='n', xlim=c(1,size), ylim=c(0,1))
    apply(net.list[[each.chr]], 1, function(x) segments(as.numeric(x['Start.bp']), as.numeric(x['LOH']), 
                                                        as.numeric(x['End.bp']), as.numeric(x['LOH']), col="red"))
  }
  dev.off()
}
