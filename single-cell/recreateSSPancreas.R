## Recreating the single-cell Pancreas visualization used
## to identify alpha and beta cells.

library(tsne)
library ("GSEABase")
library ("GSVA")
library ("GSVAdata")


####  FUNCTIONS ####
annotate <- function(hgnc){x
  require("biomaRt")
  require("org.Hs.eg.db")
  
  ens <- mapIds(org.Hs.eg.db,
                keys=hgnc,
                column="ENSEMBL",
                keytype="SYMBOL",
                multiVals="first")
  
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

getRampCol <- function(gene, colid, expr.mat, alpha.scale=20, gsea=FALSE){
  gene.idx <- grep(paste0("^", gene, "(_|$)"), rownames(expr.mat))
  #gene.col <- colorRampPalette(colors = c("grey", colid))(100)
  expr <- as.numeric(unlist(expr.mat[gene.idx,])) 
  if(gsea){
    expr[expr < 0] <- 0
  } else {
    expr[expr > quantile(expr, 0.95)] <- quantile(expr, 0.95)
  }
  col.intensity <- round(expr / max(expr), 2) / alpha.scale
  alpha(colid, col.intensity)
}

#### PRE-PROCESSING ####
sspanc <- read.table(file.path("Desktop", "tmp", "ssPancreas",
                               "GSE85241_cellsystems_dataset_4donors_updated.csv"),
                     header=TRUE, stringsAsFactors = FALSE, check.names = FALSE)
corr <- 1 - cor(sspanc, method = "pearson")
panc.tsne <- tsne(corr, perplexity=30)

ens.ids <- annotate(gsub("__.*", "", rownames(sspanc)))
sspanc.mat <- as.matrix(sspanc)
rownames(sspanc.mat) <- ens.ids
ss.zeroes <- which(apply(sspanc.mat, 1, quantile, 0.95) == 0)
sspanc.mat <- sspanc.mat[-ss.zeroes,]
sspanc.mat <- sspanc.mat[-which(is.na(rownames(sspanc.mat))),]
pancreas.scores <- gsva(sspanc.mat,
                         Pancreas.gsc,
                         method = "gsva",
                         min.sz = 2,max.sz = 500,
                         ssgsea.norm = T)
pancreas.scores[,1:20]

#### PLOTTING ####
require(scales)
plot(panc.tsne, pch=16, col=alpha("grey", 0.50))
points(panc.tsne, col=getRampCol('GCG', 'red', sspanc), pch=16)
points(panc.tsne, col=getRampCol('INS', 'blue', sspanc), pch=16)
points(panc.tsne, col=getRampCol('SST', 'green', sspanc), pch=16)
points(panc.tsne, col=getRampCol('PPY', 'yellow', sspanc), pch=16)
points(panc.tsne, col=getRampCol('KRT19', 'black', sspanc), pch=16)
points(panc.tsne, col=getRampCol('PRSS1', 'purple', sspanc), pch=16)
points(panc.tsne, col=getRampCol('MAFB', 'orange', sspanc), pch=16)

plot(panc.tsne, pch=16, col=alpha("grey", 0.50))
points(panc.tsne, col=getRampCol('Alpha', 'red', pancreas.scores, gsea = TRUE), pch=16)
points(panc.tsne, col=getRampCol('Beta', 'blue', pancreas.scores, gsea = TRUE), pch=16)
points(panc.tsne, col=getRampCol('Delta', 'green', pancreas.scores, gsea = TRUE), pch=16)
points(panc.tsne, col=getRampCol('PP', 'yellow', pancreas.scores, gsea = TRUE), pch=16)


## manually viewing descriptive stats of gsea
genes <- c('GCG', 'INS', 'SST', 'PPY', 'KRT19', 
           'PRSS1', 'MAFB')
gene.idx <- sapply(genes, function(x) grep(paste0("^", x, "_"), 
                                           rownames(sspanc)))
ord.sspanc <- rbind(sspanc[gene.idx,],
                    pancreas.scores)
ecdf.panc <- apply(ord.sspanc, 1, ecdf)
threshold <- 0.9
gene.class <- sapply(names(ecdf.panc), function(row.id){
  hi.thresh <- ecdf.panc[[row.id]](ord.sspanc[row.id, ]) > threshold
  hi.thresh[hi.thresh] <- gsub("__.*", "", row.id)
  hi.thresh[hi.thresh == 'FALSE'] <- ''
  hi.thresh
})
gene.class <- apply(gene.class[,1:4], 1, paste, collapse="")
gene.class[gene.class==''] <- NA
ord.sspanc <- rbind(ord.sspanc, gene.class)
ord.sspanc <- as.data.frame(t(ord.sspanc))
spl.sspanc <- split(ord.sspanc, f=ord.sspanc$`12`)
sspanc.summ <- lapply(spl.sspanc, summary)

sspanc.summ[[1]]
