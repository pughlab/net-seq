sample.cutoff <- 10
pow.col <- 'darkgrey'
save.env <- 0
alpha.val <- 0.15



pan_can.file <- '/Users/rquevedo/Desktop/PughLab/NET-seq/absolute_pan_can/nbt.2203-S2.txt'

# Read in absolute segment files formatted to include metadata
pan.can <- read.csv(pan_can.file,header = TRUE, sep = "\t", check.names=FALSE)
print(paste("Removing ", length(which(is.na(pan.can$ploidy))), " samples that contain no ploidy value", sep=""))
pan.can <- as.data.frame(do.call("rbind", apply(pan.can, 1, function(x) if(!is.na(x['ploidy'])) x)))
pan.can.list <- split(pan.can, pan.can$PRIMARY_DISEASE)


# Calculate the median for total_loh per disease subtype
ploidy.median <- sapply(pan.can.list, function(x) median(as.numeric(as.character(x$ploidy))))
ploidy.median[order(ploidy.median)]
pan.can.df <- do.call("rbind", pan.can.list[names(ploidy.median[order(ploidy.median)])])
pan.can.df$ploidy <- as.numeric(as.character(pan.can.df$ploidy))
pan.can.df$purity <- as.numeric(as.character(pan.can.df$purity))
pan.can.df$PRIMARY_DISEASE = factor(pan.can.df$PRIMARY_DISEASE,unique(pan.can.df$PRIMARY_DISEASE))


# # Create the fraction of LOH  to total_loh
# disease.cn.mean <- sapply(split(all.frac.loh.df,  all.frac.loh.df$disease), function(x) getLohRatio(x$cn_loh, x$tot_loh))
# disease.del.mean <- sapply(split(all.frac.loh.df,  all.frac.loh.df$disease), function(x) getLohRatio(x$del_loh, x$tot_loh))
# disease.aneu.mean <- sapply(split(all.frac.loh.df,  all.frac.loh.df$disease), function(x) getLohRatio(x$aneu_loh, x$tot_loh))
# 
# t.tot.loh.frac <- t(data.frame(cn=disease.cn.mean,
#                                del=disease.del.mean,
#                                aneu=disease.aneu.mean))

###############################
#           Visualization
###############################
box.x.val <- boxplot(ploidy~PRIMARY_DISEASE, data=pan.can.df, ylim=c(0,10))
sample.size <- box.x.val$n
# sample.pos <- which(box.x.val$names %in% sample.id)   #Label your own sample

pdf("cn_pancan.pdf")
layout(matrix(c(0,3,2,1,1,1,1,4,4,4,0), nrow=11, byrow=TRUE))
# PLOT1: Boxplot of the quantile for total LOH genomic fraction
par(mar=c(0, 4.1, 0, 4.1))
boxplot.col <- rep(pow.col, length(sample.size))
boxplot.col[which(sample.size < sample.cutoff)] <- alpha("grey", alpha.val)
#boxplot.col[sample.pos] <- 'orange'
boxplot(ploidy~PRIMARY_DISEASE, data=pan.can.df, 
        ylim=c(0,10), ylab="ploidy",
        xaxt='n', yaxt='n', las=2,
        col=boxplot.col)
axis(4, yaxp=c(0,10,10), las=2)

# PLOT 2/3: Barplot for number of samples up to 100
par(mar=c(0, 4.1, 0.5, 4.1))
x.val <- box.x.val$n 
x.val[which(x.val >= 100)] <- 100
barplot(x.val, ylim=c(0,100), yaxp=c(0,100, 2), ylab="Samples", las=1, col=boxplot.col)
abline(h=sample.cutoff, col="red")

#Barplot for number of samples above 100
par(mar=c(0.5,4.1,2.1,4.1))
x.val <- box.x.val$n 
x.val[which(x.val < 100)] <- 0
barplot(x.val, ylim=c(100,400), yaxp = c(100, 500, 1), las=1, col=pow.col)


# PLOT 4: Stacked barplot for copy-loss, copy-neutral and copy-gain distribution of total_loh
par(mar=c(5.1, 4.1, 0, 4.1))
barplot.col <- rep("black", length(sample.size))
barplot.col[which(sample.size < sample.cutoff)] <- alpha("black", alpha.val)
barplot.col[sample.pos] <- 'orange'

boxplot(purity~PRIMARY_DISEASE, data=pan.can.df, 
        yaxp = c(0, 1, 5), las=2,
        xaxt='n',
        ylab="purity",
        col=boxplot.col)
dis.labels <- as.vector(unique(pan.can.df$PRIMARY_DISEASE))
dis.labels[which(sample.size < sample.cutoff)] <- ''
axis(1, at=c(1:length(dis.labels)), labels=dis.labels, las=2, col.axis='black', cex.axis=0.7)
dis.labels <- as.vector(unique(pan.can.df$PRIMARY_DISEASE))
dis.labels[which(sample.size >= sample.cutoff)] <- ''
axis(1, at=c(1:length(dis.labels)), labels=dis.labels, las=2, col.axis=alpha('black', (alpha.val+0.2)), cex.axis=0.7)
# dis.labels <- colnames(t.tot.loh.frac)
# dis.labels[-sample.pos] <- ''
# axis(1, at=bar.x, labels=dis.labels, las=2, col.axis='orange', cex.axis=0.7)
# Add colours for each column

dev.off()



###### Creating environment
if(save.env == 1){
  pancan.cn.env <- new.env()
  pancan.cn.env$pan_can_file <- pan.can
  pancan.cn.env$pre_graph_df <- pan.can.df
  save(pancan.cn.env, file = "pancan_cn_env.Rda")
}
###### \Creating environment