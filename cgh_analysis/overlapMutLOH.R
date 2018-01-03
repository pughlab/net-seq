library(intervals)
library(reshape2)
library(ggplot2)

seq.dir <- '/Users/rquevedo/Desktop/PughLab/NET-seq/Sequenza/all'
mut.dir <- '/Users/rquevedo/Desktop/PughLab/NET-seq/Mutation_calls/mutect_haplotypecaller/filt'

gg_color_hue <- function(n) {
  hues = seq(15, 375, length = n + 1)
  hcl(h = hues, l = 65, c = 100)[1:n]
}

net.ids <- paste("NET", c("003", "009", "002", "006", "007"), "T_1", sep="-")
net.ids <- c(net.ids, paste("NET", c("001", "003", "008", "009", "002", "006", "007"), "T_2", sep="-"))

cn.mut.list <- lapply(net.ids, function(each.id){
  print(each.id)
  seg.df <- read.table(file.path(seq.dir, paste0(each.id, "-N.gz_segments.txt")), sep="\t",
                       stringsAsFactors=FALSE, check.names=FALSE, header=TRUE)
  mut.df <- read.table(file.path(mut.dir, paste0(each.id, ".snp.indels.maf")), sep="\t", 
                       header=TRUE, quote="", comment.char = "#", fill = NA, 
                       stringsAsFactors = FALSE, check.names=TRUE)
  seg.df$length <- with(seg.df, end.pos - start.pos)
  seg.df$LOH <- 0
  seg.df[which((seg.df$A == 0 | seg.df$B == 0) & seg.df$CNt > 0), 'LOH'] <- 1
  mut.list <- split(mut.df, mut.df$Chromosome)
  seg.list <- split(seg.df, seg.df$chromosome)
  
  names(mut.list) <- paste0("chr", names(mut.list))
  if(any('chrM' == names(mut.list))) mut.list <- mut.list[-grep('chrM', names(mut.list))]
  
  for(chr.id in names(mut.list)){
    mut.interval <- Intervals(with(mut.list[[chr.id]], matrix(c(Start_position, End_position), ncol=2)))
    seg.interval <- Intervals(with(seg.list[[chr.id]], matrix(c(start.pos, end.pos), ncol=2)))
    cn.stat <- do.call("rbind", lapply(interval_overlap(mut.interval, seg.interval),
                                       function(x){
                                         if(length(x) == 0) return(NA)
                                         seg.list[[chr.id]][x, c("CNt", "A", "B", "LOH"),drop=FALSE]
                                         #seg.list[[chr.id]][unlist(x), ,drop=FALSE]
                                       }))
    if(is.na(cn.stat)) cn.stat <- matrix(rep(NA, 4), nrow=1)
    colnames(cn.stat) <- c("CNt", "A", "B", "LOH")
    mut.list[[chr.id]] <- cbind(mut.list[[chr.id]], cn.stat)
  }
  mut.df <- do.call("rbind", mut.list)
  return(list("mutCN"=mut.df, 
              "segCN"=seg.df))
})

names(cn.mut.list) <- net.ids
lohhet.mut.mat <- sapply(cn.mut.list, function(x){
  x <- x[['mutCN']]
  if(any(is.na(x$LOH))) x <- x[-which(is.na(x$LOH)),]
  tbl.loh <- table(x$LOH) / nrow(x)
  if(dim(tbl.loh) == 1){
    tbl.loh <- cbind(tbl.loh, 0)
    colnames(tbl.loh) <- c(0, 1)
  }
  return(tbl.loh)
  })
rownames(lohhet.mut.mat) <- c("het", "loh")
melt.mut.lohhet <- melt(lohhet.mut.mat)
melt.mut.lohhet$metric <- 'Mut*'

seg.loh.mat <- sapply(cn.mut.list, function(x){
  x <- x[['segCN']]
  if(any(is.na(x$LOH))) x <- x[-which(is.na(x$LOH)),]
  
  x$length <- as.numeric(as.character(x$length))
  genome.size <- sum(x$length)
  tbl.loh <- matrix(c(sum(x[which(x$LOH == 0), 'length']) / genome.size,
                      sum(x[which(x$LOH == 1), 'length']) / genome.size), nrow=1)
  colnames(tbl.loh) <- c(0,1)
  return(tbl.loh)
})
rownames(seg.loh.mat) <- c("het", "loh")
melt.seg.lohhet <- melt(seg.loh.mat)
melt.seg.lohhet$metric <- 'Seg**'

melt.seg.mut <- rbind(melt.mut.lohhet, melt.seg.lohhet)
colnames(melt.seg.mut) <- c("het_loh", "sid", "value", "metric")
melt.seg.mut <- melt.seg.mut[order(melt.seg.mut$sid),]

loh.frac.df <- as.data.frame(cbind(seg.loh.mat['loh',], lohhet.mut.mat['loh',]))
colnames(loh.frac.df) <- c("cn_seg", "mut_load")

pdf("/Users/rquevedo/Desktop/PughLab/NET-seq/Mutation_calls/mutect_haplotypecaller/mutation_loh_skew.pdf",
    height = 3)
ggplot(melt.seg.mut, aes(x = metric, y = value, fill = het_loh)) + 
  geom_bar(stat = 'identity', position = 'stack') + facet_grid(~ sid) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  labs(x="",y="Fraction") 
dev.off()

pdf("/Users/rquevedo/Desktop/PughLab/NET-seq/Mutation_calls/mutect_haplotypecaller/mutation_loh_skew.corr.pdf", 
    width=5, height=5)
plot(loh.frac.df, col=gg_color_hue(2)[2], pch=16, 
     ylim=c(0,1), xlim=c(0,1), las=1, asp=1,
     xlab="Genomic LOH fraction", ylab="Mutational LOH fraction")
abline(lm(data=loh.frac.df, mut_load ~ cn_seg))

dev.off()

cor.test(loh.frac.df[,1], loh.frac.df[,2], method='pearson')$p.val
cor.test(loh.frac.df[,1], loh.frac.df[,2], method='pearson')$estimate




