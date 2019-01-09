### Analysis of single-cell Pancreas data:  Project the 
### SRR and NET-seq samples on to the ssPancreas data.

## LIBRARIES
library(AnnotationDbi)
library(TxDb.Hsapiens.UCSC.hg19.knownGene)
library(org.Hs.eg.db)
library(biomaRt)

## VARIABLES
pdir <- '/mnt/work1/users/pughlab/projects/NET-SEQ/rna_seq_external/clustering'

#### FUNCTIONS #### 
# annotates a vector of HGNC ids using two methods
annotate <- function(hgnc, input='hgnc'){
  require("biomaRt")
  require("org.Hs.eg.db")
  
  if('input'=='hgnc'){
    ens <- mapIds(org.Hs.eg.db,
                  keys=hgnc,
                  column="ENSEMBL",
                  keytype="SYMBOL",
                  multiVals="first")
  } else {
    ens <- mapIds(org.Hs.eg.db,
                  keys=hgnc,
                  column="SYMBOL",
                  keytype="ENSEMBL",
                  multiVals="first")
  }

  if(any(is.na(ens))){
    mart <- useMart("ENSEMBL_MART_ENSEMBL")
    mart <- useDataset("hsapiens_gene_ensembl", mart)
    na.ens <- getBM(mart=mart, 
                    attributes=c("external_gene_name", "hgnc_symbol", "ensembl_gene_id"),
                    filter="external_gene_name",
                    values=names(ens[which(is.na(ens))]),
                    uniqueRows=TRUE)
    ids <- sapply(na.ens$external_gene_name, function(x) {
      grep(paste0("^", x, "$"), names(ens), ignore.case=TRUE)
    })
    ens[unlist(ids)] <- rep(na.ens$ensembl_gene_id,
                            sapply(ids, length))
  }
  
  ens
}


#### MAIN #### 
# Reading in the single cell pancreas data
setwd(pdir)
sspanc <- read.table(file.path("input", "GSE85241_cellsystems_dataset_4donors_updated.csv"),
                     header=TRUE, stringsAsFactors = FALSE, check.names = FALSE)
anno.df <- data.frame(hugo=gsub("__.*", "", rownames(sspanc)),
                      stringsAsFactors = FALSE) 
anno.df$ensembl <- annotate(anno.df$hugo)
sspanc$ensembl <- anno.df$ensembl

# Reading and combining annotation files
load(file.path(pdir, "input", 
               list.files(file.path(pdir, "input"), 
                          pattern = "^SRR.*annotation.*Rdata")[1]))
global.anno.df$ensembl <- gsub("\\..*$", "", global.anno.df$gene_id)
global.anno.df <- merge(global.anno.df, anno.df, by='ensembl', all.x=TRUE)

# Reading and combining the tpm mats with single cell
tpm.mats <- lapply(list.files(file.path(pdir, "input"), 
                  pattern = "^(SRR|NET).*tpm.*Rdata"), 
       function(each.tpm){
         load(file.path(pdir, "input", each.tpm))
         tpm.mat <- as.data.frame(tpm.mat)
         tpm.mat$ensembl <- gsub("\\..*$", "", rownames(tpm.mat))
         tpm.mat
})
tpm.mats[[length(tpm.mats) + 1]] <- sspanc
tpm.mat <- Reduce(function(x, y) merge(x, y, by = "ensembl", all.x = T), 
                  tpm.mats)
tpm.mat <- as.matrix(tpm.mat[,-1])
rownames(tpm.mat) <- global.anno.df$gene_id

# Marking SS ids
ss.ids <- colnames(tpm.mats[[length(tpm.mats)]])
ss.ids <- ss.ids[-length(ss.ids)]

## ssGSEA with single cell gene sigs
load(file.path("ref", "ssPancreas.Rdata"))
tpm.tmp <- tpm.mat
rownames(tpm.tmp) <- gsub("\\..*", "", rownames(tpm.tmp))
ss.zeroes <- which(apply(tpm.tmp, 1, quantile, 0.95) == 0)
tpm.tmp <- tpm.tmp[-ss.zeroes,]
#tpm.tmp <- tpm.tmp[-which(is.na(rownames(tpm.tmp))),]
pancreas.scores <- gsva(tpm.tmp,
                        Pancreas.gsc,
                        method = "gsva",
                        min.sz = 2,max.sz = 500,
                        ssgsea.norm = T)

# PCA to nominate variable genes
tpm.mat[is.na(tpm.mat)] <- 0
t.tpm.mat <- t(tpm.mat)
ss.tpm.mat <- t.tpm.mat[(rownames(t.tpm.mat) %in% ss.ids),]
ss.zeroes <- which(apply(ss.tpm.mat, 2, quantile, 0.95) == 0)
t.tpm.mat <- t.tpm.mat[,-ss.zeroes]

scaled.tpm.mat <- scale(t.tpm.mat)
pca.tmp <- prcomp(scaled.tpm.mat)
highvar.pc <- which((cumsum(pca.tmp$sdev) / sum(pca.tmp$sdev)) <= 0.95)
comp <- data.frame(pca.tmp$x[,highvar.pc])
tau.cor <- cor(comp, method = "kendall",  use = "pairwise")
rho.cor <- cor(t(comp), method = "spearman",  use = "pairwise")
pearson.cor <- cor(t(comp), method="pearson", use="pairwise") 

pc.hc <- hclust(dist(comp))
pc.hc2 <- hclust(as.dist(1-pearson.cor))

## Get IDs of single cells
g.ids <- data.frame(id=rownames(t.tpm.mat), class="ssCell", 
                    stringsAsFactors = FALSE)
g.ids[grep("^NET-", g.ids$id), 'class'] <- 'NETSEQ'
g.ids[grep("^SRR", g.ids$id), 'class'] <- 'SRR'

hgnc <- annotate(gsub("\\..*$", "", colnames(t.tpm.mat)), 
                 input='ensembl')
marker.genes <- c('Alpha'="GCG", 'Beta'="INS", 
                  'Delta'='SST', 'PP'='PPY')
cell.type.idx <- match(marker.genes, hgnc)
ss.idx <- which(g.ids$class == 'ssCell')

marker.mat <- t.tpm.mat[ss.idx, cell.type.idx]
#z.scores <- apply(marker.mat, 2, scale)
z.scores <- apply(marker.mat, 2, function(x) ecdf(x)(x))
max.z <- apply(z.scores, 1, which.max)
low.z <- which(apply(z.scores,1, max) < 0.75)
max.marker <- names(marker.genes)[max.z]
max.marker[low.z] <- 'ssCell'
g.ids[ss.idx, 'class'] <- max.marker

library(ggplot2)
library(dendextend)
library(RColorBrewer)
n <- 60
qual_col_pals = brewer.pal.info[brewer.pal.info$category == 'qual',]
col_vector = unlist(mapply(brewer.pal, qual_col_pals$maxcolors, rownames(qual_col_pals)))
col.mapping <- data.frame(col=col_vector[1:length(unique(g.ids$class))],
                          class=unique(g.ids$class))
g.ids$cols <- col.mapping[match(g.ids$class, col.mapping$class), 'col']

pc.hc2$labels <- as.character(g.ids$class)
change.idx <- grep("(SRR|NETSEQ)", pc.hc2$labels, invert = TRUE)
pc.hc2$labels[change.idx] <- "."

dend.pc <- pc.hc2 %>% as.dendrogram %>%
  set("branches_lwd", 0.5) %>%
  set("labels_colors") %>%  
  set("leaves_pch", 19) %>% 
  set("leaves_col", g.ids$cols[pc.hc2$order])

ggd1 <- as.ggdend(dend.pc)

pdf("dend.ss.pdf")
leg <-  unique(g.ids[,c('class', 'cols')])
plot(0, type='n', ylim=c(1,150), xlim=c(1,50))
legend("topleft", legend=leg$class, fill=leg$cols, cex=0.8)

ggplot(ggd1, labels = TRUE) + 
  scale_y_reverse(expand = c(0.2, 0)) +
  coord_polar(theta="x")
dev.off()