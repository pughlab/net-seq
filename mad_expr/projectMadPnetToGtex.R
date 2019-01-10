## LIBRARIES

## VARIABLES
pdir <- '/mnt/work1/users/pughlab/projects/NET-SEQ/rna_seq_external/clustering'

## MAIN
setwd(pdir)

# Loads in all TPM matrices and records their IDs/class
loaded.data <- lapply(list.files(file.path(pdir, "input"), 
                  pattern = "tpm.*Rdata"), function(tpm.id){
  uid <- gsub(".tpm.*", "", tpm.id)
  print(paste0("Loading in ", uid, "..."))
  anno.id <- paste0(uid, ".annotation.Rdata")
  load(file.path("input", anno.id))
  load(file.path("input", tpm.id))
  
  tpm.mat <- round(tpm.mat, 2)
  id.class <- data.frame(class=rep(uid, ncol(tpm.mat)),
                         id=colnames(tpm.mat))
  
  list("tpm"=tpm.mat, "id"=id.class)
})

# Aggregate and GC the TPM data and annotations
g.tpm.mat <- do.call("cbind", lapply(loaded.data, function(x) x[['tpm']]))
g.ids <- do.call("rbind", lapply(loaded.data, function(x) x[['id']]))


g.tpm.mat[is.na(g.tpm.mat)] <- 0
t.tpm.mat <- t(g.tpm.mat)
zeroes <- which(apply(t.tpm.mat, 2, quantile, 0.95) == 0)
t.tpm.mat <- t.tpm.mat[,-zeroes]


scaled.tpm.mat <- scale(t.tpm.mat)
pca.tmp <- prcomp(scaled.tpm.mat)
highvar.pc <- which((cumsum(pca.tmp$sdev) / sum(pca.tmp$sdev)) <= 0.95)
comp <- data.frame(pca.tmp$x[,highvar.pc])
rho.cor <- cor(t(comp), method="spearman", use="pairwise") 
pearson.cor <- cor(t(comp), method="pearson", use="pairwise") 
pc.hc2 <- hclust(as.dist(1-pearson.cor))
pc.hc <- hclust(dist(comp))

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
change.idx <- grep("(Liver|Pancreas|SRR|NETSEQ)", pc.hc2$labels, invert = TRUE)
pc.hc2$labels[change.idx] <- "."

dend.pc <- pc.hc2 %>% as.dendrogram %>%
  set("branches_lwd", 0.5) %>%
  set("labels_colors") %>%  
  set("leaves_pch", 19) %>% 
  set("leaves_col", g.ids$cols[pc.hc2$order])

ggd1 <- as.ggdend(dend.pc)

pdf("dend.pdf")
leg <-  unique(g.ids[,c(1,3)])
plot(0, type='n', ylim=c(1,150), xlim=c(1,50))
legend("topleft", legend=leg$class, fill=leg$cols, cex=0.8)

ggplot(ggd1, labels = TRUE) + 
  scale_y_reverse(expand = c(0.2, 0)) +
  coord_polar(theta="x")
dev.off()