library("GEOquery")
library(dplyr)
require(GSVA)
require(GSEABase)
require(msigdbr)
require(TxDb.Hsapiens.UCSC.hg19.knownGene)
library(org.Hs.eg.db)
require(BSgenome.Hsapiens.UCSC.hg19)
require(matrixStats)
require(scales)

###################
#### FUNCTIONS ####
cleanMetadata <- function(gsm, meta.file=NULL, group=NULL){
  if(group=='sadanadam'){
    ## Mine the IDs by mixing GSM data with an external worksheet
    metadata <- read.table(meta.file, sep="\t", 
                           header=TRUE, stringsAsFactors = FALSE, check.names=FALSE)
    
    metadata <- as.data.frame(metadata) 
    colnames(metadata)[c(1,6)] <- c("sample", "mad")
    
    ## Format metadata into MAD groups
    cl <- sapply(gsm, function(x) strsplit(x@header$title, "-")[[1]][1])
    cl <- t(sapply(seq_along(cl), function(i) c(cl[i], names(cl)[i])))
    cl <- as.data.frame(cl); colnames(cl) <- c("sample", "gsm")
    metadata <- merge(metadata[,c(1,6)], cl, by='sample')
    
    ## Relabel into M1, D, and A annotation, similar used in Chan
    metadata$mad <- sapply(strsplit(metadata$mad, ","), function(i){
      mad.bool <- c(any(grepl("DAXX", i)), 
                    any(grepl("ATRX|ARTX", i, perl=TRUE)), 
                    any(grepl("MEN1", i)))
      wt.stat <- if(!any(mad.bool) & length(i) != 0) 'WT' else NA
      mad <- c("D", "A", "M1")
      if(!is.na(wt.stat)) wt.stat else paste(mad[mad.bool], collapse="")
    })
    metadata <- metadata[-which(metadata$mad == ''),]
    
    cl.spl <- split(metadata, f=metadata$mad)
    cl.spl[['MAD']] <- do.call("rbind", cl.spl[names(cl.spl) != 'WT'])
  } else if(group=='chan'){
    ## Mine the ID's from the GSM data
    classification <- sapply(gsm, function(x) strsplit(x@header$title, "_")[[1]][1])
    cl.spl <- split(classification, f=classification)
    mad.ids <- unlist(cl.spl[names(cl.spl) != 'WT'])
    names(mad.ids) <- gsub("^.*\\.", "", names(mad.ids))
    cl.spl[['MAD']] <- mad.ids
  } else {
    stop("Group must be either 'sadanadam' or 'chan'")
  }
  
  cl.spl
}

grpExprMat <- function(gse, cl.spl, group=NULL){
  # Clean up multiple probe mappings
  .cleanExprMat <- function(e.gse, id.df){
    # Find genes with multiple probes and reduce them using sample mean
    mult.genes <- names(which(table(id.df[,2])>1))
    comb.idx <- sapply(mult.genes, function(i) {
      grep(paste0("^", i, "$"), x=id.df[,2], ignore.case=TRUE)
    })
    na.idx <- which(is.na(id.df[,2]))
    
    keep.idx <- c(1:nrow(id.df))[-c(na.idx, unlist(comb.idx))]
    KEEP.e.gse <- e.gse[keep.idx,]
    rownames(KEEP.e.gse) <- id.df[keep.idx,2]
    COMB.e.gse <- t(sapply(comb.idx, function(i) colMeans(e.gse[i,])))
    new.e.gse <- rbind(KEEP.e.gse, COMB.e.gse)
    
    new.e.gse
  }
  
  e.gse <- exprs(gse)
  if(group=='sadanadam'){
    rownames(e.gse) <- featureData(gse)@data$ORF
  } else if(group == 'chan'){
    require("hgu133plus2.db")
    ids <- as.list(hgu133plus2SYMBOL[featureNames(gse)])
    locs <- as.list(hgu133plus2CHRLOC[featureNames(gse)])
    id.df <- t(sapply(seq_along(ids), function(i) c(names(ids)[i], ids[[i]])))
    e.gse <- .cleanExprMat(e.gse, id.df)
  } else {
    stop("group must be either 'sadanadam' or 'chan'")
  }
  
  expr.spl <- lapply(cl.spl, function(grp){
    if(group == 'chan') g.id = names(grp) else g.id = grp$gsm
    idx <- which(colnames(e.gse) %in% g.id)
    
    e.gse[,idx]
  })
  
  expr.spl
}

compareGoi2Cen <- function(goi, e.gse, expr.spl, gen.plot=FALSE){
  goi.idx <- grep(goi, rownames(e.gse))
  goi.exprs <- lapply(expr.spl, function(i) i[goi.idx, ,drop=FALSE])
  
  # Compare GOI to itself across groups
  goi.summ <- sapply(goi.exprs, function(i) summary(as.numeric(i)))
  goi.ttest <- sapply(goi.exprs, function(a) {
    round(sapply(goi.exprs, function(b) t.test(a,b)$statistic),3)
  })
  
  ## Check the correlation of expression between GOI and Centromeric proteins
  cenp.id <- "CENP|MAD2|BUB|KIF2C|AURK"
  cenp.ids <- rownames(e.gse)[grep(cenp.id, rownames(e.gse), perl=TRUE)]
  cor.mats <- lapply(cenp.ids, function(cenp){
    cenp.idx <- grep(cenp, rownames(e.gse))
    cenp.exprs <- lapply(expr.spl, function(i) i[cenp.idx[1], ,drop=TRUE])
    
    sapply(names(goi.exprs), function(cls){
      if(ncol(goi.exprs[[cls]]) > 2){
        ct <- cor.test(goi.exprs[[cls]][1,],
                       cenp.exprs[[cls]], method='pearson')
        round(c(ct$p.value, ct$estimate),3)
      } else { c(NA, NA) }
      
    })
  })
  names(cor.mats) <- cenp.ids
  
  cor.pval <- t(sapply(cor.mats, function(i) i[1,]))
  cor.estimate <- t(sapply(cor.mats, function(i) i[2,]))
  
  ## Visualization
  if(gen.plot) {
    boxplot(goi.exprs, ylim=c(0,max(sapply(goi.exprs, max))))
    boxplot(t(sapply(cor.mats, function(i) i[2,])))
  }
  
  list("n"=sapply(expr.spl, ncol),
       "goi"=list("ttest"=goi.ttest,
                  "summ"=goi.summ),
       "cor"=list("pval"=cor.pval,
                  "estimate"=cor.estimate))
}

genGsc <- function(keyterms, categ='C5', subcat='BP'){
  m = msigdbr(species = "Homo sapiens",
              category=categ, 
              subcategory = subcat)
  
  ## Find all genesets related to the keyterms
  goterms <- sapply(keyterms, function(go){
    unique(m[grep(go, m$gs_name, ignore.case=TRUE),]$gs_name)
  })
  gsoi <- lapply(unlist(goterms), function(go){
    m[grep(go, m$gs_name, ignore.case=TRUE),]$human_gene_symbol  
  })
  names(gsoi) <- unlist(goterms)
  
  ## Generate a geneset collection
  centromere.gsc <- GeneSetCollection( 
    lapply(names(gsoi), function(x) {
      GeneSet(setName=x,
              geneIds=unique(sort(gsoi[[x]])),
              geneIdType =SymbolIdentifier(),
              collectionType = ExpressionSetCollection())
    }))
  
  centromere.gsc
}

grpAndTestSSGSEA <- function(expr.spl, gsc){
  ## Calculate ssGSEA scores
  ssgsea.scores <- lapply(expr.spl, function(expression){
    centromere.scores <- gsva(expression,
                              gsc,
                              method = "ssgsea",
                              min.sz = 2,max.sz = 500,
                              ssgsea.norm = T)
    centromere.scores
  })
  
  ## Test for group differences
  t.mat <- sapply(seq_along(rownames(ssgsea.scores[[1]])), function(r.idx){
    t.val <- t.test(ssgsea.scores[['MAD']][r.idx,], ssgsea.scores[['WT']][r.idx,])
    c("stat"=t.val$statistic, "p"=t.val$p.value)
  })
  colnames(t.mat) <- rownames(ssgsea.scores[[1]])
  
  ## Stats of each geneset 
  genes <- sapply(centromere.gsc@.Data, function(i) geneIds(i))
  names(genes) <- sapply(centromere.gsc@.Data, setName)
  
  list("score.test"=t(t.mat),
       "geneset"=list("n"=sapply(genes, length),
                      "genes"=genes))
}

mapExprsGeneToLoci <- function(e.gse){
  require(TxDb.Hsapiens.UCSC.hg19.knownGene)
  require(org.Hs.eg.db)
  txdb = TxDb.Hsapiens.UCSC.hg19.knownGene # abbreviate
  gene.gr <- genes(txdb)
  gene.gr$symbol <- select(org.Hs.eg.db, 
                           keys=gene.gr$gene_id, 
                           keytype="ENTREZID", 
                           columns="SYMBOL")$SYMBOL
  gene.gr <- as.data.frame(gene.gr)
  
  expr.genes <- rownames(e.gse)
  anno.expr.genes <- gene.gr[match(expr.genes, gene.gr$symbol),
                             c("seqnames", "start", "end", "symbol")]
  
  order.idx <- with(anno.expr.genes, order(seqnames, start))
  list("order"=order.idx,
       "anno"=anno.expr.genes[order.idx,])
}

plotSSGSEA <- function(ssgsea, grp.col){
  orig.go <- rownames(ssgsea)
  rownames(ssgsea) <- c("CEN", paste0("KIN", 1:4))
  
  idx <- barplot(ssgsea[,1], col=alpha(grp.col, 0.5), 
                 xlim=c(-5, 5), border = NA, las=1,
                 xlab='t-statistic', ylab='', horiz = TRUE)
  
  ystar.idx <- idx[which(ssgsea[,2] <= 0.15 & ssgsea[,2] > 0.05)]
  xstar.idx <- ssgsea[which(ssgsea[,2] <= 0.15 & ssgsea[,2] > 0.05),1]
  ytwostar.idx <- idx[which(ssgsea[,2] <= 0.05)]
  xtwostar.idx <- ssgsea[which(ssgsea[,2] <= 0.05),1]
  
  text(x = xstar.idx, y=ystar.idx, labels = "*", 
       pos = if(xstar.idx < 0) 2 else 4, col="black", srt=90)
  text(x = xtwostar.idx, y=ytwostar.idx, labels = "**", 
       pos = if(xtwostar.idx < 0) 2 else 4, col="black", srt=90)
  
  data.frame("orig"=orig.go,
             "new"=rownames(ssgsea))
}

plotExprZscore <- function(alt.exprs, ref.exprs, 
                           anno.genes, grp.col, ...){
  alt.exprs <- alt.exprs[anno.genes[['order']],]
  ref.exprs <- ref.exprs[anno.genes[['order']],]
  
  ## Calculate z-score for all MAD versus WT samples
  .calcZ <- function(x){ apply(x, 1, scale) }
  all.exprs <- cbind(alt.exprs, ref.exprs)
  z.all.exprs <- .calcZ(all.exprs)
  z.alt.exprs.mu <- colMeans(z.all.exprs[1:ncol(alt.exprs),])
  z.alt.exprs.sd <- colSds(z.all.exprs[1:ncol(alt.exprs),])
  
  z.alt <- data.frame("low"=z.alt.exprs.mu - 1*z.alt.exprs.sd,
                      "mu"=z.alt.exprs.mu,
                      "hi"=z.alt.exprs.mu + 1*z.alt.exprs.sd)
  
  .getChrLength <- function(genome = "BSgenome.Hsapiens.UCSC.hg19"){
    g <- getBSgenome(genome, masked=FALSE)
    seqlengths(g)[1:24]
  }
  chr.len <- .getChrLength()
  chr.len <- rbind(chr.len, c(0, cumsum(as.numeric(chr.len)))[-25])
  
  anno.chr <- split(anno.genes[['anno']], f=anno.genes[['anno']]$seqnames)[1:24]
  z.chr <- split(z.alt, f=anno.genes[['anno']]$seqnames)[1:24]
  
  plot(0, type='n', xlim=c(0, sum(chr.len[1,])),
       xaxt='n', ylab='Z', xlab='', las=2, ...)
  abline(v=chr.len[2,], lty=2)
  axis(side = 1, at = (chr.len[2,] + (chr.len[1,]/2)), 
       labels = colnames(chr.len), 
       las=2, cex.axis=0.7, tick = FALSE)
  
  sapply(names(anno.chr)[1:24], function(chr){
    chr.spacer <- chr.len[2,grep(paste0("^", chr, "$"), colnames(chr.len))]
    gene.x <- with(anno.chr[[chr]], chr.spacer + start)
    z.chr[[chr]]$index <- gene.x
    loessHiZ50 <- loess(hi ~ index, data=z.chr[[chr]], span=0.50) # 50% smoothing span
    loessLoZ50 <- loess(low ~ index, data=z.chr[[chr]], span=0.50) # 50% smoothing span
    loessMuZ50 <- loess(mu ~ index, data=z.chr[[chr]], span=0.50) # 50% smoothing span
    
    polygon(x = c(loessLoZ50$x, 
                  rev(loessHiZ50$x)), 
            y = c(loessLoZ50$fitted,
                  rev(loessHiZ50$fitted)),
            border=NA,
            col=alpha(grp.col, 0.5))
    lines(x=loessMuZ50$x, loessMuZ50$fitted, col=grp.col, lwd=2)
    
  })
  
}
###################
#### VARIABLES ####
pdir <- '/mnt/work1/users/pughlab/projects/NET-SEQ/hgu133a2_array_chan/expr_analysis'


###############################
#### Sadanadam et al, 2015 ####
setwd("/mnt/work1/users/pughlab/projects/NET-SEQ/external_data/Sadanandam_exprArray")
grp.col <- '#d95f02'

## Load Datasets
gse=getGEO(filename="GSE73338_series_matrix.txt.gz")
soft=getGEO(filename="GSE73338_family.soft.gz")
gsm <- soft@gsms

meta.file <- "144026_2_supp_3166825_nvh38t.txt"
cl.spl <- cleanMetadata(gsm, meta.file, 'sadanadam')
expr.spl <- grpExprMat(gse, cl.spl, group='sadanadam')
e.gse <- exprs(gse)

compareGoi2Cen('DAXX', e.gse, expr.spl, gen.plot=FALSE)
as.matrix(sort(t(sapply(cor.mats, function(i) i[1,]))[,'MAD']))
as.matrix(sort(t(sapply(cor.mats, function(i) i[2,]))[,'MAD']))

centromere.gsc <- genGsc(c("centromere", "kinetochore"))
ssgsea.stats <- grpAndTestSSGSEA(expr.spl[c("MAD", "WT")], centromere.gsc)

pdf(file.path(pdir, "sad_ssgseaMAD.pdf"), width = 6, height = 4)
plotSSGSEA(ssgsea.stats[['score.test']], grp.col)
dev.off()

pdf(file.path(pdir, "sad_exprMAD.pdf"), width = 7, height = 5)
anno.genes <- mapExprsGeneToLoci(e.gse)
plotExprZscore(expr.spl[['MAD']], expr.spl[['WT']],
               anno.genes, grp.col, ylim=c(-1, 1))
dev.off()


##########################
#### Chan et al, 2018 ####
setwd("/mnt/work1/users/pughlab/projects/NET-SEQ/external_data/Chan-hgu133a2")
grp.col <- '#7570b3'

gse=getGEO(filename="GSE117851_series_matrix.txt.gz")
soft=getGEO(filename="GSE117851_family.soft.gz")
gsm <- soft@gsms

cl.spl <- cleanMetadata(gsm, group='chan')
expr.spl <- grpExprMat(gse, cl.spl, group='chan')
e.gse <- exprs(gse)

compareGoi2Cen('DAXX', e.gse, expr.spl, gen.plot=FALSE)
as.matrix(sort(t(sapply(cor.mats, function(i) i[1,]))[,'MAD']))
as.matrix(sort(t(sapply(cor.mats, function(i) i[2,]))[,'MAD']))

centromere.gsc <- genGsc(c("centromere", "kinetochore"))
ssgsea.stats <- grpAndTestSSGSEA(expr.spl[c("MAD", "WT")], centromere.gsc)

pdf(file.path(pdir, "chan_ssgseaMAD.pdf"), width = 6, height = 4)
plotSSGSEA(ssgsea.stats[['score.test']], grp.col)
dev.off()

pdf(file.path(pdir, "chan_exprMAD.pdf"), width = 7, height = 5)
anno.genes <- mapExprsGeneToLoci(e.gse)
plotExprZscore(expr.spl[['MAD']], expr.spl[['WT']],
               anno.genes, grp.col, ylim=c(-1, 1))
dev.off()