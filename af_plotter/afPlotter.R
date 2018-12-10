.libPaths(c(.libPaths(), "/mnt/work1/users/pughlab/bin/r_lib/library/3.1"))
###############################################################
#
#                 Allelic Fraction Plotter v2.0
#   Purpose: Bins the SNPs in WGS data into X-size bins
#     takes the median and CBS across those segments to 
#     estimate purity based on LOH segments
###############################################################
library(VariantAnnotation)
library(scales)
library(data.table)
library(intervals)
library(copynumber)
library(DNAcopy)

##############################
##### Setting up the workspace
# Set the variables
tumor.norm <- FALSE   # If tumor and normal are given, and bins are computed, fits distributions and tests for difference
purity.limit <- FALSE  # If ydiscovering the lower purity limit is of interest
purity.bins <- TRUE
purity.points <-FALSE
rm.x <- FALSE
global.bin.size <- 50000  # genomic bin size
global.chrom <- c(1:22) # chromosome ID
assembly <- 'GRCh37'  # GRCh37 or GRCh38
purity.levels <- c(0.25, 0.50, 0.70, 0.80)
cn.stat <- 2     # assumed copy-neutral ploidy for LOH segments, default = 2
project.id <- 'chan-adm'

# Processing variables
snvs <- 0   #whether mutect SNVs are given
min.depth <- 15   # minimum depth to consider a snp

# Load the RData and functions
load('~/git/net-seq/af_plotter/data/afPlotter.ref.Rdata') # Contains: clinvar.df, refgene.df, chrom.df
source('~/git/net-seq/af_plotter/src/binProcess.R')
source('~/git/net-seq/af_plotter/src/snvProcess.R')
source('~/git/net-seq/af_plotter/src/snpProcess.R')
source('~/git/net-seq/af_plotter/src/visualizationAf.R')
source('~/git/net-seq/af_plotter/src/compareSeg.R')
source('~/git/net-seq/af_plotter/src/misc.R')

#Load in the SNPs
if(project.id %in% 'netseq'){
  # RQ NET-Seq Project:
  mutect.onco.maf <- '~/git/net-seq/af_plotter/data/snvs'
  exomeseq.vcf <- '~/git/net-seq/af_plotter/data/snps_wes'
  rnaseq.vcf <- '~/git/net-seq/af_plotter/data/snps_rna'
  snvs <- 1   #whether mutect SNVs are given
  
  rna.net.001 <- "NET-001-T_4_STAR_Aligned_rg_added_sorted_dups_marked.bam_Haplotype_caller_out_VariantFiltration.vcf.gz"
  rna.net.002 <- 'NET-002-T_4_STAR_Aligned_rg_added_sorted_dups_marked.bam_Haplotype_caller_out_VariantFiltration.vcf.gz'
  rna.net.003 <- "NET-003-T_4_STAR_Aligned_rg_added_sorted_dups_marked.bam_Haplotype_caller_out_VariantFiltration.vcf.gz"
  rna.net.006 <- "NET-006-T_4_STAR_Aligned_rg_added_sorted_dups_marked.bam_Haplotype_caller_out_VariantFiltration.vcf.gz"
  rna.net.007 <- "NET-007-T_4_STAR_Aligned_rg_added_sorted_dups_marked.bam_Haplotype_caller_out_VariantFiltration.vcf.gz"
  rna.net.008 <- "NET-008-T_4_STAR_Aligned_rg_added_sorted_dups_marked.bam_Haplotype_caller_out_VariantFiltration.vcf.gz"
  rna.net.009 <- "NET-009-T_4_STAR_Aligned_rg_added_sorted_dups_marked.bam_Haplotype_caller_out_VariantFiltration.vcf.gz"
  
  dna.net.001 <- "NET-001-T_2.snps.indelx.vcf"
  dna.net.002 <- "NET-002-T_2.snps.indelx.vcf"
  dna.net.003 <- "NET-003-T_2.snps.indelx.vcf"
  dna.net.003a <- 'NET-003-T_1.snps.indelx.vcf'
  dna.net.006 <- "NET-006-T_2.snps.indelx.vcf"
  dna.net.007 <- "NET-007-T_2.snps.indelx.vcf"
  dna.net.008 <- "NET-008-T_2.snps.indelx.vcf"
  dna.net.009 <- "NET-009-T_2.snps.indelx.vcf"
  dna.net.009a <- 'NET-009-T_1.snps.indelx.vcf'
  
  group1 <- matrix(c(dn1=dna.net.001, rn1=rna.net.001,     # PNET Samples
                     dn3=dna.net.003, rn3=rna.net.003,
                     dn3a=dna.net.003a, rn3=rna.net.003,
                     dn8=dna.net.008, rn8=rna.net.008,
                     dn9=dna.net.009, rn9=rna.net.009,
                     dn9a=dna.net.009a, rn9=rna.net.009), byrow=T, ncol=2)
  colnames(group1) <- c("dna", "rna")
  group2 <-  matrix(c(dn2=dna.net.002, rn2=rna.net.002,   # GINET Samples
                      dn6=dna.net.006, rn6=rna.net.006,
                      dn7=dna.net.007, rn7=rna.net.007), byrow=T, ncol=2)
  colnames(group2) <- c("dna", "rna")
  
  mutect.net.001 <- 'NET-001-T_2.snp.indels.maf'
  mutect.net.003 <- 'NET-003-T_2.snp.indels.maf'
  mutect.net.003a <- 'NET-003-T_1.snp.indels.maf'
  mutect.net.008 <- 'NET-008-T_2.snp.indels.maf'
  mutect.net.009 <- 'NET-009-T_2.snp.indels.maf'
  mutect.net.009a <- 'NET-009-T_1.snp.indels.maf'
  group1 <- cbind(group1, c(mutect.net.001, mutect.net.003, mutect.net.003a, 
                            mutect.net.008, mutect.net.009, mutect.net.009a))
  colnames(group1) <- c("dna", "rna", "snv")
  
  mutect.net.002 <- 'NET-002-T_2-N.snp.indels.maf'
  mutect.net.006 <- 'NET-006-T_2-N.snp.indels.maf'
  mutect.net.007 <- 'NET-007-T_2-N.snp.indels.maf'
  group2 <- cbind(group2, c(mutect.net.002, mutect.net.006, mutect.net.007))
  colnames(group2) <- c("dna", "rna", "snv")
  
  dir.create("~/Desktop/netseq/net_disc/output/", recursive = TRUE, showWarnings = FALSE)
  output.dir <- '~/Desktop/netseq/net_disc/output/'
} else if(project.id %in% 'chan-adm'){
  PDIR = '/mnt/work1/users/pughlab/projects/NET-SEQ/rna_seq_external/af_plots'
  source(file.path(PDIR, "scripts", "plotAF.R"))
  
  rnaseq.vcf <- file.path(PDIR, "input")
  exomeseq.vcf <- file.path(PDIR, "input")
  output.dir <- file.path(PDIR, "output")
  snvs <- 0   #whether mutect SNVs are given
  
  vcf.files <- list.files(rnaseq.vcf, "vcf$")
  group1 <- data.frame("dna"=vcf.files,
                       "rna"=vcf.files)
}


# Generate preliminary reference dataframes and genomic bins
clinvar.df <- formatClinvarDf(clinvar.df, assembly)
chrom.bins.list <- generateGenomeBins(ucsc.chrom.df ,global.bin.size)


##############################
##### Processing the SNPs
sample.match.m <- group1
rownames(sample.match.m) <- gsub("\\..+", "", sample.match.m[,'dna'])
rnasegs <- list()
for(each.sample in 1:dim(sample.match.m)[1]){
  # Load in the SNPs from VCF and separate based on chromosome
  dna.af.list <- generateChrList(sample.match.m[each.sample,'dna'], exomeseq.vcf)
  rna.af.list <- generateChrList(sample.match.m[each.sample,'rna'], rnaseq.vcf)
  
  # Remove X Chromosomes:
  if(rm.x){
    dna.af.list[['chrX']] <- NULL
    rna.af.list[['chrX']] <- NULL
  }
  
  # Annotate based on mutect SNVs
  if(dim(sample.match.m)[2] > 2){
    dna.af.snv.list <- annotateSnvs(sample.match.m[each.sample, 'snv'], mutect.onco.maf, dna.af.list)
    rna.af.snv.list <- annotateSnvs(sample.match.m[each.sample, 'snv'], mutect.onco.maf, rna.af.list)
    
    dna.af.snv.list <- removeHomSnps('dna', dna.af.snv.list, rna.af.snv.list)
    rna.af.snv.list <- removeHomSnps('rna', dna.af.snv.list, rna.af.snv.list)
  } else {
    dna.af.snv.list <- lapply(dna.af.list, function(x) {
      x$Variant_Classification <- NA
      x$Hugo_Symbol <- NA
      return(x)
    })
    rna.af.snv.list <- lapply(rna.af.list, function(x) {
      x$Variant_Classification <- NA
      x$Hugo_Symbol <- NA
      return(x)
    })
  }
  
  #Overlap and maintain pathogenic snps/snvs
  dna.af.list.ol <- overlapAf('dna', dna.af.snv.list, rna.af.snv.list)
  rna.af.list.ol <- overlapAf('rna', dna.af.snv.list, rna.af.snv.list)
  dna.af.list.ol <- dna.af.snv.list
  rna.af.list.ol <- rna.af.snv.list
  
  #Bin the SNPs into x-sized genomic bins to accomodate for noise
  dna.af.list.bin <- binSnps(dna.af.list.ol, chrom.bins.list, rm.hom=FALSE, shallow.wgs=FALSE)
  #rna.af.list.bin <- dna.af.list.bin
  rna.af.list.bin <- binSnps(rna.af.list.ol, chrom.bins.list, rm.hom=FALSE, shallow.wgs=FALSE)

  dir.create(file.path(output.dir, rownames(sample.match.m)[each.sample]), recursive=TRUE)
  setwd(file.path(output.dir, rownames(sample.match.m)[each.sample]))  
  # Retrieve the Allelic fraction determined by CBS segmentation for each SNP
  dna.rna.list <- plotAfExomeRna(paste(rownames(sample.match.m)[each.sample], ".ol", sep=""), 
                                 exome.list=dna.af.list.ol, rna.list=rna.af.list.ol,
                                 ebin.list=dna.af.list.bin, rbin.list=rna.af.list.bin,
                                 min.depth=15, filt.region=TRUE, 
                                 plot.bins=FALSE, plot.points=TRUE, bin.size=global.bin.size)
  rnasegs[[rownames(sample.match.m)[each.sample]]] <- dna.rna.list[['dna.seg']]
}

save(rnasegs, file=file.path(output.dir, "rnasegs.RData"))

rnasegs.df <- lapply(names(rnasegs), function(id) {
  do.call("rbind", lapply(rnasegs[[id]], function(chr.m) {
    chr.m$ID <- id
    chr.m
  }))
})
write.table(do.call("rbind", rnasegs.df), file=file.path(output.dir, "rnasegs.seg"),
            sep="\t", quote = FALSE , col.names = TRUE, row.names = FALSE)