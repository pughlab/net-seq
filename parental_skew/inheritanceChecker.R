##### This code find the correlations between samples from a vcf input file. It is made for Sample check #####
### The vcf file : The vcf file is created by merging all samples' vcfs from unifiedgenotyper using vcftools## the merged vcf has all the genotype data after the 10th column###
library(scales)
rm(list =ls())
###########
#arguments##
###########
sample_dir='/mnt/work1/users/pughlab/projects/NET-SEQ/exome_seq/qc/parental_GenotypeChecker/input'
sample_mega_file='merge.vcf'

###########
#Functions#
###########
#### Genotype reader
genotype <- function(b) {
  c = b[[1]][1]
  if (c == "." | c == "0/0"){genotype = "0/0"}
  else {genotype = c}
  return(genotype)}

######
#Main##
#######
setwd(file.path(sample_dir, "allChr"))
### Read the Sample genotype table and assign columns###
Mega_VCF = read.csv(sample_mega_file, sep = "\t", comment.char = "#")
print(paste("All sample names are: ", 
            colnames(Mega_VCF[,c(10:ncol(Mega_VCF))]), sep=" "))

rownames(Mega_VCF) = paste0(Mega_VCF$CHROM,":", Mega_VCF$POS, 
                            Mega_VCF$REF, ">", Mega_VCF$ALT)
sample_sub = Mega_VCF[,10:ncol(Mega_VCF)]

# Parse haplotype caller into the first column (i.e. 0/0, or 1/1, or 0/1).  Stores as dataframe
sample_geno <- apply(sample_sub, 2, function(x) strsplit(as.vector(x), split = ":"))
sample_count <- lapply(sample_geno, function(i) sapply(i, genotype))
sample_count.data <- as.data.frame(do.call("cbind", sample_count))
sample_count.data[] <- lapply(sample_count.data, function(x) as.character(x))

# Trios index
p.idx <- c("parental"=colnames(sample_count.data)[4],
           "normal"=colnames(sample_count.data)[1],
           "tumor"=colnames(sample_count.data)[3])

# for(x in c(2:10)){
#   sample_count.data[,p.idx['tumor']] <- gsub(x, "1", sample_count.data[,p.idx['tumor']])
# }

# Because the tumor stat is obscured by purity issues, I estimate the allelic fraction
# based on ref to alt depth.
tumor.df <- do.call("rbind", sample_geno[[p.idx['tumor']]])
t.ref.alt.depth <- strsplit(tumor.df[,5], split=',')  # index 5 is the ref, alt depth
t.homstat <- sapply(t.ref.alt.depth, function(x){
  af <- tryCatch({
    round(as.integer(x[1]) / (as.integer(x[1]) + as.integer(x[2])),3)
  }, error=function(e){NA})
  
  depth <- tryCatch({
    as.integer(x[1]) + as.integer(x[2])
  }, error=function(e){NA})
  
  stat <- tryCatch({ 
    stat <- 'het'
    # If the ref to alt depth is below 0.2 or above 0.8, i term is 'ref' or 'alt' accordingly
    if(af <= 0.2 | af >= 0.8){
      if(as.integer(x[1]) > as.integer(x[2])) stat <- 'ref' else stat <- 'alt'
    } 
    stat
  }, error=function(e){ stat <- 'het'})
  return(c("af"=af, "stat"=stat, "depth"=depth))
})
t.homstat <- t(t.homstat)
# Hard set all 'refs' to 0/0 and all 'alts' to 1/1 since they are being called "het" 
# because of purity contamination
sample_count.data[which(t.homstat[,2] == 'ref'), p.idx['tumor']] <- '0/0'
sample_count.data[which(t.homstat[,2] == 'alt'), p.idx['tumor']] <- '1/1'

# Identify all homozygous snps in parent and heterozygous in normal
homP.idx <- grep("0/0|1/1", sample_count.data[,p.idx['parental']])
hetN.idx <- grep("0/1", sample_count.data[,p.idx['normal']])
# Find homozygous parent to heterozygous in child
hom.to.het.idx <- intersect(homP.idx, hetN.idx)
# Identify whether these heterozygous SNPs are now matching parental or not in the tumour
match.idx <- sample_count.data[hom.to.het.idx, p.idx['tumor']] == sample_count.data[hom.to.het.idx, p.idx['parental']]
overlap.match <- table(match.idx)[2]/length(hom.to.het.idx) # HOM match / All Hom to Het SNPs
# number of heterozygous SNPs in tumor sample 
het.tumor.cnt <- table(sample_count.data[hom.to.het.idx,p.idx['tumor']] == '0/1')[2]
mat <- matrix(c(round(overlap.match,3), 
                length(hom.to.het.idx),
                het.tumor.cnt), ncol=3)
## num.snps: Number of homozygous SNPs in parent that are heterozygous in offspring
## overlap.frac: Fraction of num.snps that are homozygous in the tumour
## num.het: Number of num.snps that are still heterozygous in tumour
colnames(mat) <- c("overlap.frac", "num.snps", "num.het")
print(mat)

pdf("AF.refAlt.pdf")
p.col <- rep(alpha("red", 0.50), length(hom.to.het.idx))  # non-parental red
p.col[match.idx] <- alpha('blue', 0.50) # match parental - blue
het.idx <- which(sample_count.data[hom.to.het.idx,p.idx['tumor']] == '0/1')
p.col[het.idx] <- alpha('black', 0.5) # heterozygous SNP - black

plot(x=Mega_VCF[hom.to.het.idx, 'POS'], 
     y=as.numeric(t.homstat[hom.to.het.idx,1]), 
     col=p.col, ylim=c(-0.1, 1), xlab="Position", ylab="AF")
points(x=Mega_VCF[hom.to.het.idx, 'POS'], y=rep(-0.05, length(hom.to.het.idx)), 
       pch=15, col=p.col)

dev.off()

write.table(mat, file = "sampleMatch.tsv",
            row.names = TRUE, col.names = NA , sep = "\t")

