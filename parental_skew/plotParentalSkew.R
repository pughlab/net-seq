library(scales)
## Fig Xa
## Visualize the parental skew of tumour SNPs per chromosome arm
## Output is from parentalidChecker.R run on each chromosome
## The output from each chromosome was aggregated together
#
# rm allChr_sampleMatch.tsv
# 
# paste -d "\t" <(echo "Chr") <(head -1 chr1.sampleMatch.tsv) > allChr_sampleMatch.tsv
# 
# paste -d "\t" <(ls -1 | grep "\.sampleMatch.tsv" | sed s/".sampleM.*"/""/) \
# <(for i in $(ls -1 | grep "\.sampleMatch.tsv"); do tail -n+2 $i; done) \
# >> allChr_sampleMatch.tsv
#
## num.snps: Number of homozygous SNPs in parent that are heterozygous in offspring
## overlap.frac: Fraction of num.snps that are homozygous in the tumour
## num.het: Number of num.snps that are still heterozygous in tumour

PDIR <- '/mnt/work1/users/pughlab/projects/NET-SEQ/exome_seq/qc/parental_GenotypeChecker/input/allChr'
file <- 'allChr_sampleMatch.tsv'

geno <- read.table(file.path(PDIR, file), sep="\t",
                   header=TRUE, stringsAsFactors = FALSE)
chrs <- paste0("chr", c(1:22, "X", "Y"))
geno <- geno[order(factor(geno$Chr, levels=chrs)),]
geno$het.frac <- round(with(geno, num.het / num.snps),3)
geno$paternal.frac <- round(with(geno, 1 - (overlap.frac + het.frac)),3)
colnames(geno) <- c('Chr', 'X', 'maternal.frac', 'num.snps', 
                    'num.het', 'het.frac', 'paternal.frac')
rownames(geno) <- geno$Chr

chr.ord <- rev(rownames(geno))

outdir <- gsub("input", "output", PDIR)
dir.create(outdir, recursive = TRUE, showWarnings = FALSE)
pdf(file.path(outdir, "parental_skew.pdf"), width = 5)
split.screen(matrix(c(0, 0.6, 0, 1,
                      0.6, 1.0, 0, 1), byrow = TRUE, ncol=4))
screen(1); par(mar=c(5.1, 4.1, 4.1, 0.4))
barplot(t(geno[,c("maternal.frac", "het.frac", "paternal.frac")])[,chr.ord],
        col=c("maroon1", alpha("grey", 0.5), "navyblue"),
        horiz = TRUE, las=1, xlab='Fraction of SNPs', xaxt='n', 
        border = NA)
axis(side = 1, at = c(0, 0.5, 1.0), labels = c(0, 0.5, 1.0))
screen(2); par(mar=c(5.1, 0.3, 4.1, 2.1))
barplot(t(geno[chr.ord, 'num.snps', drop=FALSE]), horiz=TRUE, 
        las=1, yaxt='n', xlab='SNPs (x1000)', xaxt='n',
        border = NA)
axis(side = 1, at = c(0, 1000, 2000), labels = c(0, 1, 2))
close.screen(all.screens = TRUE); dev.off()

write.table(x = geno, file = file.path(outdir, "parental_skew.tsv"),
            sep="\t", quote = FALSE, row.names = FALSE, col.names = TRUE)