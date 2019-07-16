pdf("~/Desktop/misseg_Loh.pdf", height = 3)

r <- biserial.cor(missegregateChr()[1:22], lohChr()[1:22])
loh.cols <- lohChr()
loh.cols[loh.cols=='black'] <- "#33a02c"
loh.cols[loh.cols=='white'] <- "#fb9a99"

## Plot the misseg. chromosomes as a function of LOH chromosomes
set.seed(1)
pdf("~/Desktop/lohMisseg.pdf")
plot(missegregateChr()[1:22],
     as.integer(factor(lohChr()[1:22])) + runif(n = 22, min=-0.3, max=0.3),
     col=alpha(loh.cols, 0.6), pch=16,
     yaxt='n', xlab="Misseg. Fraction", ylab="", ylim=c(0, 5), axes=FALSE)
axis(side = 1, at = seq(0, 1, by=0.2), labels=seq(0, 1, by=0.2))
axis(side=2, at=c(1,2), labels=c("LOH", "Het"), tick = FALSE, las=1)
text(x=1, y=2.5, labels=paste0("r = ", round(r,3), " (biserial)"), pos=2)
dev.off()


## Plot Missegregate chr and LOH chromosomes in a barplot
loh.bar <- abs(as.integer(factor(lohChr())) - 2)
names(loh.bar) <- names(lohChr())
pdf("~/Desktop/misseg_Loh.pdf", height = 3)
split.screen(c(2,1))
screen(1); par(mar=c(0.5, 6.1, 4.1, 2.1))
barplot(missegregateChr()[1:22], col="black", las=2, 
        xaxt='n', yaxt='n', ylab="Misseg.\nfraction")
axis(side = 2, at=c(0, 0.5, 1.0), labels = c(0, 0.5, 1.0), las=2)
sapply(c(1:4), axis, at=c(-100, 100))
abline(v=x[-length(x)] + (diff(x) / 2))

screen(2); par(mar=c(5.1, 6.1, 0.5, 2.1))
x <- barplot(loh.bar[1:22], col="black", las=2, 
             yaxt='n', ylab="LOH\nstatus")
sapply(c(1:4), axis, at=c(-100, 100))
abline(v=x[-length(x)] + (diff(x) / 2))
close.screen(all.screens=TRUE)
dev.off()