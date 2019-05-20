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

require(beeswarm)
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

grpExprMat <- function(gse, cl.spl, group=NULL, return.tmat=FALSE){
  # Clean up multiple probe mappings
  .cleanExprMat <- function(e.gse, gene.ids){
    # Find genes with multiple probes and reduce them using sample mean
    mult.genes <- names(which(table(gene.ids)>1))
    comb.idx <- sapply(mult.genes, function(i) {
      grep(paste0("^", i, "$"), x=gene.ids, ignore.case=TRUE)
    })
    na.idx <- which(is.na(gene.ids))
    
    keep.idx <- c(1:nrow(e.gse))[-c(na.idx, unlist(comb.idx))]
    KEEP.e.gse <- e.gse[keep.idx,]
    rownames(KEEP.e.gse) <- gene.ids[keep.idx]
    COMB.e.gse <- t(sapply(comb.idx, function(i) colMeans(e.gse[i,])))
    new.e.gse <- rbind(KEEP.e.gse, COMB.e.gse)
    
    new.e.gse
  }
  
  e.gse <- exprs(gse)
  if(group=='sadanadam'){
    rownames(e.gse) <- featureData(gse)@data$ORF
    id.df <- rownames(e.gse)
  } else if(group == 'chan'){
    require("hgu133plus2.db")
    ids <- as.list(hgu133plus2SYMBOL[featureNames(gse)])
    id.df <- t(sapply(seq_along(ids), function(i) c(names(ids)[i], ids[[i]])))
    id.df <- id.df[,2]
  } else {
    stop("group must be either 'sadanadam' or 'chan'")
  }
  
  e.gse <- .cleanExprMat(e.gse, id.df)
  
  if(return.tmat){
    expr.spl <- e.gse
  } else {
    expr.spl <- lapply(cl.spl, function(grp){
      if(group == 'chan') g.id = names(grp) else g.id = grp$gsm
      idx <- which(colnames(e.gse) %in% g.id)
      
      e.gse[,idx]
    })
  }
  
  expr.spl
}

compareGoi2Cen <- function(goi, e.gse, expr.spl, core.kin, gen.plot=FALSE){
  goi.idx <- sapply(goi, grep, x=rownames(e.gse))
  goi.exprs <- lapply(expr.spl, function(i) i[goi.idx, ,drop=FALSE])
  
  ## Compare GOI to itself across groups
  goi.ttest <- lapply(goi.exprs, function(a) {
    lapply(goi.exprs, function(b) {
      m <- sapply(1:nrow(a), function(gene.a){
        sapply(1:nrow(b), function(gene.b){
          t.test(a[gene.a,],b[gene.b,])$statistic
        })
      })
      colnames(m) <- rownames(a); rownames(m) <- rownames(b)
      m
    })
  })
  
  ## Check the correlation of expression between GOI and Centromeric proteins
  cenp.id <- core.kin[,1]
  cenp.ids <- rownames(e.gse)[which(rownames(e.gse) %in% cenp.id)]
  #cenp.id <- "CENP|MAD2|BUB|KIF2C|AURK"
  #cenp.ids <- rownames(e.gse)[grep(cenp.id, rownames(e.gse), perl=TRUE)]
  cor.mats <- lapply(cenp.ids, function(cenp){
    cenp.idx <- match(cenp, rownames(e.gse))
    cenp.exprs <- lapply(expr.spl, function(i) i[cenp.idx[1], ,drop=TRUE])
    
    sapply(names(goi.exprs), function(cls){
      if(ncol(goi.exprs[[cls]]) > 2){
        sapply(rownames(goi.exprs[[cls]]), function(g){
          ct <- cor.test(goi.exprs[[cls]][g,],
                         cenp.exprs[[cls]], method='pearson')
          #round(c(ct$p.value, ct$estimate),3)
          round(ct$estimate,3)
        })
      } else { NULL }
    })
  })
  names(cor.mats) <- cenp.ids
  
  ## Visualization
  if(gen.plot) {
    if(length(goi) > 1) warning("This was not designed to work with more than 1 GOI")
    boxplot(goi.exprs, ylim=c(0,max(sapply(goi.exprs, max))))
    boxplot(t(sapply(cor.mats, function(i) i[2,])))
  }
  
  list("n"=sapply(expr.spl, ncol),
       "goi"=goi.ttest,
       "cor"=cor.mats)
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

plotBeeswarmCor <- function(cor.mats, goi, group='core-kinetochore', core.kin){
  mad <- as.data.frame(sapply(goi, function(g, grp){
    sapply(cor.mats[['cor']], function(i) i[paste0(g, ".cor"), grp])
  }, grp='MAD'))
  
  wt <- as.data.frame(sapply(goi, function(g, grp){
    sapply(cor.mats[['cor']], function(i) i[paste0(g, ".cor"), grp])
  }, grp='WT'))
  
  colnames(core.kin) <- c("Symbol", "ID", "group")
  .annoGeneCor <- function(ds, mad.stat){
    ds$Symbol <- rownames(ds)
    ds$mad <- mad.stat
    merge.ds <- merge(ds, core.kin, by="Symbol", all.x=TRUE)
    merge.ds$color = 'black'
    merge.ds[grep("CENP", merge.ds$Symbol),'color'] <- 'red'
    merge.ds
  }
  merged.cor <- rbind(.annoGeneCor(mad, 'MAD'),
                      .annoGeneCor(wt, 'WT'))
  ck.merged.cor <- merged.cor[which(merged.cor$group == group),]
  
  split.screen(c(1,3))
  screen(1); par(mar=c(5.1, 4, 4.1, 0.1))
  beeswarm(DAXX ~ mad, data=ck.merged.cor,breaks=NA, 
           main='DAXX', ylab='Correlation', xlab='',
           pwcol=ck.merged.cor$color, pch=16, ylim=c(-0.6, 0.6))
  screen(2); par(mar=c(5.1, 4, 4.1, 0.1))
  beeswarm(ATRX ~ mad, data=ck.merged.cor,breaks=NA, 
           main='ATRX', ylab='', xlab='', yaxt='n',
           pwcol=ck.merged.cor$color, pch=16, ylim=c(-0.6, 0.6))
  axis(side = 2, at=seq(-1, 1, by=0.2), labels = rep("", 11))
  screen(3); par(mar=c(5.1, 4, 4.1, 0.1))
  beeswarm(MEN1 ~ mad, data=ck.merged.cor,breaks=NA, 
           main='MEN1', ylab='', xlab='',  yaxt='n',
           pwcol=ck.merged.cor$color, pch=16, ylim=c(-0.6, 0.6))
  axis(side = 2, at=seq(-1, 1, by=0.2), labels = rep("", 11))
  close.screen(all.screens=TRUE)
  
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
core.file <- '/mnt/work1/users/pughlab/projects/NET-SEQ/external_data/thiru_core-kinetochore/supp_E14-03-0837_mc-E14-03-0837-s01.txt'
core.kin <- read.table(core.file, sep="\t", header=TRUE, 
                       stringsAsFactors = FALSE, check.names = FALSE)

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
e.gse <- grpExprMat(gse, cl.spl, group='sadanadam', return.tmat = TRUE)

## Generate Correlation beeswarm plots highlighting CENP genes
goi <- c("DAXX", "MEN1", "ATRX")
cor.mats <- compareGoi2Cen(goi, e.gse, expr.spl[c("MAD", "WT")],
                           core.kin, gen.plot=FALSE)
pdf(file.path(pdir, "sad_beeswarmcor.pdf"), width = 6, height = 7)
plotBeeswarmCor(cor.mats, goi, group='core-kinetochore', core.kin)
dev.off()

## Generate ssGSEA plots showing disruption of centromeric process
centromere.gsc <- genGsc(c("centromere", "kinetochore"))
ssgsea.stats <- grpAndTestSSGSEA(expr.spl[c("MAD", "WT")], centromere.gsc)

pdf(file.path(pdir, "sad_ssgseaMAD.pdf"), width = 6, height = 4)
plotSSGSEA(ssgsea.stats[['score.test']], grp.col)
dev.off()

## Generate z-score expression plots showing copy-number
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
e.gse <- grpExprMat(gse, cl.spl, group='chan', return.tmat = TRUE)

## Generate Correlation beeswarm plots highlighting CENP genes
goi <- c("DAXX", "MEN1", "ATRX")
cor.mats <- compareGoi2Cen(goi, e.gse, expr.spl[c("MAD", "WT")],
                           core.kin, gen.plot=FALSE)
pdf(file.path(pdir, "chan_beeswarmcor.pdf"), width = 6, height = 7)
plotBeeswarmCor(cor.mats, goi, group='core-kinetochore', core.kin)
dev.off()

## Generate ssGSEA plots showing disruption of centromeric process
centromere.gsc <- genGsc(c("centromere", "kinetochore"))
ssgsea.stats <- grpAndTestSSGSEA(expr.spl[c("MAD", "WT")], centromere.gsc)

pdf(file.path(pdir, "chan_ssgseaMAD.pdf"), width = 6, height = 4)
plotSSGSEA(ssgsea.stats[['score.test']], grp.col)
dev.off()

pdf(file.path(pdir, "ssgseaGeneSets.pdf"), width = 10)
plot(0, type='n', xlim=c(0,10), ylim=c(0,10), axes=FALSE, xlab='', ylab='')
text(x = 2, y=seq(1,10,by=2), labels=c("CEN", paste0("KIN", 1:4)), pos=2)
text(x = 10, y=seq(1,10,by=2), labels=names(ssgsea.stats$geneset$n), pos=2)
dev.off()

## Generate z-score expression plots showing copy-number
pdf(file.path(pdir, "chan_exprMAD.pdf"), width = 7, height = 5)
anno.genes <- mapExprsGeneToLoci(e.gse)
plotExprZscore(expr.spl[['MAD']], expr.spl[['WT']],
               anno.genes, grp.col, ylim=c(-1, 1))
dev.off()