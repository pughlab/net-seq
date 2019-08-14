library(scales)

telomere <- read.table("~/Onedrive - UHN/PughLab/NET-seq/telomere/allTelbamLengths.csv", 
                       sep=",", header=TRUE, stringsAsFactors = FALSE)
head(telomere)
telomere$Sample <- gsub("_02_.*", "", telomere$Sample)
pnet <- paste0("NET-2-", c('029', '005', '030', '031',
                           '034', '025', '003', '013', 
                           '017', '002', '011', '015', '032'))
ginet <- paste0("NET-2-", c('037', '006', '024', '027',
                            '004', '009', '021', '033',
                            '016', '007'))
net <- c(pnet, ginet)
telomere <- telomere[which(telomere$Sample %in% net),]
telo.spl <- split(telomere, f=telomere$Sample %in% pnet)
names(telo.spl) <- c('GINET', 'PNET')

pdf("~/Onedrive/PughLab/NET-seq/telomere/telomere-length.pdf", width=4)
{
boxplot(lapply(rev(telo.spl), function(x) log10(x$Length)),
        yaxt='n', col=alpha("grey", 0.5), ylim=c(3, 6),
        ylab="Telomere length (kbp)")
axis(side = 2, at = log10(c(1000, 5000, 10000, 
                            50000, 100000, 500000)), 
     labels = c(1, 5, 10, 50, 100, 500), las=2)
points(rep(2, length(ginet)), log10(telo.spl[['GINET']]$Length), 
       pch=16, col=alpha("black", 0.6))
points(rep(1, length(pnet)), log10(telo.spl[['PNET']]$Length), 
       pch=16, col=alpha("black", 0.6))
net29idx <- which(telo.spl[['PNET']] == 'NET-2-029')
points(1, log10(telo.spl[['PNET']][net29idx,'Length']), 
       pch=16, col=alpha("red", 0.8))

pval <- t.test(telo.spl[['PNET']]$Length,
               telo.spl[['GINET']]$Length, 
               alternative = "greater")$p.val
f.pval <- var.test(telo.spl[['PNET']]$Length,
                   telo.spl[['GINET']]$Length, 
                   alternative = "greater")$p.val
text(x = 1.5, y = 5.7, labels = paste0("p = ",
                                       round(pval, 4), 
                                       " (one-sided t-test)"),
     cex=0.7, pos = 3)
segments(x0 = 1,y0 = 5.7,x1 = 2,y1 = 5.7)
segments(x0 = 1,y0 = 5.6,x1 = 1,y1 = 5.7)
segments(x0 = 2,y0 = 5.6,x1 = 2,y1 = 5.7)
}
dev.off()

sapply(colnames(xy[,which(xy['Mean',] > 0.8)]), function(x){
   summary(pnet[['estq']][,x])
})