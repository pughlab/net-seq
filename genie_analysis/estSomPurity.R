library(intervals)
library(limma)
#Dataset obtained and results generated on April 27-2017
dir.create("~/Desktop/netseq/genie/plots/", recursive = TRUE, showWarnings = FALSE)
setwd("~/Desktop/netseq/genie")
setwd("~/git/net-seq/genie_analysis")
source("~/git/net-seq/genie_analysis/src/processing.R")
source("~/git/net-seq/genie_analysis/src/afCalculation.R")


####################
######  Functions
{
  getSnpAf <- function(mut.df){
    mut.df$alt.af <- with(mut.df, t_alt_count / t_depth)
    mut.df$ref.af <- with(mut.df, t_ref_count / t_depth)
    mut.split <- split(mut.df, mut.df$tumor_sample_barcode)
    mut.snp.split <- lapply(mut.split, function(x) {
      dbsnp.idx <- grep("^rs", x$dbsnp_rs)
      if(length(dbsnp.idx) > 0){
        dbsnp.df <- x[dbsnp.idx,,drop=FALSE]
        novel.df <- x[-dbsnp.idx,,drop=FALSE]
      } else {
        dbsnp.df <- NA
        novel.df <- x
      }
      return(list("dbsnp"=dbsnp.df,
                  "novel"=novel.df))
    })
  }
  
  getPurityRange <- function(x){
    min.val <- apply(x, 1, function(i) min(i, na.rm=TRUE))
    max.val <- apply(x, 1, function(i) max(i, na.rm=TRUE))
    
    return(c("min.purity" = max(min.val),
             "max.purity" = min(max.val)))
  }
}

####################
######  Preprocess Data
{
  base.colnames <- c("Tumor_Sample_Barcode", "Hugo_Symbol", "Chromosome", "Start_Position", "End_Position",
                     "Variant_Classification", "Variant_Type", "Reference_Allele", "Tumor_Seq_Allele2", "dbSNP_RS")
  base.colnames <- tolower(base.colnames)
  
  # MuTect GENIE format
  mut.df <- read.table("~/git/net-seq/genie_analysis/genie_data/GENIE_pnet_DAM_mutations.txt", header=TRUE, sep="\t",
                       check.names=FALSE, stringsAsFactors = FALSE)
  cn.df <- read.table("genie_data/GENIE_pnet_DAM_cna.seg", header=TRUE, sep="\t",
                      check.names=FALSE, stringsAsFactors = FALSE)
  #HaplotypeCaller format
  otb.mut.df <- read.table("~/git/net-seq/genie_analysis/otb_data/OTB_NETseq_haplotypecaller.maf", header=TRUE, sep="\t",
                           check.names=FALSE, stringsAsFactors = FALSE,  quote="")
  otb.cn.df <- read.table("~/git/net-seq/genie_analysis/otb_data/OTB_NETseq_varscan_cna.seg", header=TRUE, sep="\t",
                          check.names=FALSE, stringsAsFactors = FALSE)
  ploidy.models <- read.table("~/git/net-seq/genie_analysis/data/ploidy_models.tsv", sep="\t",
                              check.names=FALSE, header=TRUE)
  ploidy.models <- as.data.frame(t(ploidy.models))
  cn.df$length <- with(cn.df, loc.end - loc.start) / 1000000
  
  colnames(otb.mut.df) <- tolower(colnames(otb.mut.df))
  otb.mut.df$tumor_sample_barcode <- otb.mut.df$matched_norm_sample_barcode
  otb.mut.df <- otb.mut.df[,c(base.colnames,'allelic_depth', 'read_depth')] 
  otb.mut.df$t_alt_count <- as.integer(gsub("^.*,", "", otb.mut.df$allelic_depth))
  otb.mut.df$t_ref_count <- as.integer(gsub(",.*$", "", otb.mut.df$allelic_depth))
  otb.mut.df$t_depth <- otb.mut.df$read_depth
  otb.mut.snp.split <- getSnpAf(otb.mut.df)
  colnames(otb.cn.df) <- c("ID", "chrom", "loc.start", "loc.end", "num.mark", 
                           "normal_depth", "tumor_depth", "seg.mean", "gc_content")
  otb.cn.df$seg.mean <- 2^(otb.cn.df$seg.mean)
  otb.cn.df$length <- with(otb.cn.df, loc.end - loc.start)
  otb.cn.df$chrom <- gsub("chr", "", otb.cn.df$chrom)
  otb.cn.split <- split(otb.cn.df, otb.cn.df$ID)
  
  colnames(mut.df) <- tolower(colnames(mut.df))
  mut.df <- mut.df[,c(base.colnames, "exon_number", 
                      "t_depth", "t_ref_count", "t_alt_count",
                      "sift", "polyphen", "exac_af" )]
  mut.snp.split <- getSnpAf(mut.df)
 
  
  cn.split <- split(cn.df, cn.df$ID)
  
  
  mut.snp.split <- append(mut.snp.split, otb.mut.snp.split)
  cn.split <- append(cn.split, otb.cn.split)
  
  mut.snp.cnOV.split <- mut.snp.split[names(cn.split)]
}



####################
######  Visualization
# Visualize the CN .seg overlaid with mutational data and their AF for MSK data (tumor-normal match)
# Blue = SNVs matching rs IDs
# Red = Novel SNVs
{
  pdf("~/Desktop/netseq/genie/plots/seg_af_plots.pdf", width=16, height = 5)
  id <- 1
  net.identifiers <- data.frame()
  for(each.sample in names(cn.split)){
    print(paste0("Sample: ", each.sample))
    each.cn.df <- cn.split[[each.sample]]
    each.cn.df$seg.mean <- with(each.cn.df, seg.mean - median(with(each.cn.df, rep(seg.mean, length))))
    each.dbsnp.df <- mut.snp.cnOV.split[[each.sample]]$dbsnp
    each.novel.df <- mut.snp.cnOV.split[[each.sample]]$novel
    split.screen(matrix(c(0, 1, 0, 0.7,
                          0, 1, 0.7, 1), byrow=TRUE, ncol=4))
    screen(2)
    par(mar=c(0, 4.1, 4.1, 2.1))
    plot(0, type='n', ylab='', xlab='', yaxt='n', xaxt='n', axes=FALSE,
         ylim=c(0,10), xlim=c(0,100))
    text(x=40, y=5, labels = paste(each.sample, ":", paste0("GENIE_Pnet_", id)))
    net.identifiers <- rbind(net.identifiers, 
                             data.frame(each.sample, paste0("GENIE_Pnet_", id)))
    id <- id + 1
    
    screen(1)
    screen.idx <- split.screen(c(1, 23 + 2))
    
    screen(3)
    par(mar=c(5.1, 3.15, 0, 0)) 
    plot(0, type='n', xlim=c(0,0), ylim=c(-1.5, 1.5), 
         xlab='', ylab='', axes=FALSE)
    axis(side = 2, at=seq(-1.5, 1.5, by=0.5), labels=seq(-1.5, 1.5, by=0.5), 
         cex=0.7, las=2)
    mtext("seg.mean", side=2, line=2, cex.lab=1,las=3)
    
    chrom.idx <- c(1:22, "X")
    for(each.chrom in chrom.idx){
      chr.cn.df <- each.cn.df[grep(paste0("^", each.chrom, "$"), each.cn.df$chrom),]
      chr.dbsnp.df <- tryCatch({each.dbsnp.df[grep(paste0("^", each.chrom, "$"), each.dbsnp.df$chromosome),]},
                               error=function(e){data.frame()})
      chr.novel.df <- tryCatch({each.novel.df[grep(paste0("^", each.chrom, "$"), each.novel.df$chromosome),]},
                               error=function(e){data.frame()})
      
      print(dim(chr.novel.df))
      
      screen(match(each.chrom, chrom.idx) + 1 + 2)
      par(mar=c(5.1, 0, 0, 0))
      
      plot(0, type='n', xlim=c(0, max(chr.cn.df$loc.end)), ylim=c(-1.5, 1.5), 
           xaxt='n', yaxt='n', xlab=each.chrom)
      rect(xleft = chr.cn.df$loc.start, ybottom = chr.cn.df$seg.mean - 0.01,
           xright = chr.cn.df$loc.end, ytop = chr.cn.df$seg.mean + 0.01, col="black")
      if(nrow(chr.dbsnp.df) > 0) plotPoints(chr.dbsnp.df, chr.cn.df, "blue", 'dbsnp')
      if(nrow(chr.novel.df) > 0) plotPoints(chr.novel.df, chr.cn.df, "red", 'novel', increase.interval.size = TRUE)
    }
    close.screen(all.screens=TRUE)
  }
  dev.off()
}
net.identifiers <- data.frame("genie"=names(mut.snp.cnOV.split), 
                              "newid"=paste0("GENIE_Pnet_", c(1:length(mut.snp.split))))


# Cycles through each sample and extracts the allelic fraction and identifiers for each
mut.snp.cnOV.split <- mut.snp.split
all.pur.est <- lapply(net.identifiers$newid, function(net.id){
  print(net.id)
  af.df <- mut.snp.cnOV.split[[net.identifiers[match(net.id, net.identifiers$newid),1]]]
  
  af.list <- list()
  #Extract AF and identifiers for "dbsnp" variants
  try(expr=for(each.symbol in af.df$dbsnp$hugo_symbol){
    af.list[[paste0(each.symbol, "_snp")]] <- af.df$dbsnp[which(each.symbol ==  af.df$dbsnp$hugo_symbol),
                                                          'alt.af']
  }, silent=TRUE)  
  
  #Extract AF and identifiers for "novel" variants
  try(expr=for(each.symbol in af.df$novel$hugo_symbol){
    af.list[[each.symbol]] <- af.df$novel[which(each.symbol ==  af.df$novel$hugo_symbol),
                                                         'alt.af']
    }, silent=TRUE)
  
  mult.idx <- which(sapply(af.list, length) > 1)
  if(length(mult.idx) > 0) {
    for(each.idx in mult.idx){
      af.list[[each.idx]] <- mean(af.list[[each.idx]])
    }
  }
  
  # Runs both dbsnp and novel variants to estimate Purity across all mutations and all possible CN profiles
  genDf(af.list, id=net.id)
})
# Collapses and summarizes the data for excel
names(all.pur.est) <- net.identifiers$newid
x <- do.call("rbind", all.pur.est)
write.table(x, file="~/Desktop/netseq/genie/otb_purity_estimate.tsv", sep="\t", col.names =TRUE, quote=FALSE)
write.table(t(ploidy.models), file="~/Desktop/netseq/genie/ploidy_models.tsv", sep="\t", col.names =TRUE, quote=FALSE)


# sample.cn.profile <- data.frame("ID"=x$ID,
#                                 "Gene"=x$Gene,
#                                 "CN"=c(1, 1, 3, 1, 1),
#                                 "Chr"=c(2, 2, 2, 2, 2))
cn.profile <- read.table("~/git/net-seq/genie_analysis/data/Genie_purity_estimations2.txt", sep="\t", header=TRUE,
                        stringsAsFactors = FALSE, check.names = FALSE)
cn.profile <- split(cn.profile, cn.profile$ID)
best.fit.models <- list()
all.pur.est <- all.pur.est[c(9:37)]
names(all.pur.est) <- paste0("GENIE_Pnet_", c(1:27))
for(net.id in names(all.pur.est)){
  print(net.id)
  sample.pur <- all.pur.est[[net.id]]
  sample.cn.profile <- cn.profile[[net.id]]
  sample.pur <- sample.pur[,-c(1,2,3)]
  
  # Establish the range of acceptable purities and remove any violators
  pur.range <- getPurityRange(sample.pur)
  #sample.pur[sample.pur > pur.range['max.purity'] | sample.pur < pur.range['min.purity']] <- NA
  if(all(is.na(sample.pur))) next
  
  #Index genes
  gene.idx <- sapply(sample.cn.profile$Gene, function(x) match(x, rownames(sample.pur)))
  names(gene.idx) <- sample.cn.profile$Gene
  
  # Generate all possible purity models irrespective of copy-number/ploidy
  all.models <- list()
  for(each.gene in sample.cn.profile$Gene){
    # Preprocess purity vals and filter out autosomal/allosomal purities accordingly
    pur.vals <- sample.pur[gene.idx[each.gene],]
    normal.ploidy <- which(ploidy.models$Pn != sample.cn.profile[match(each.gene, sample.cn.profile$Gene), 'Chr'])
    pur.vals[normal.ploidy] <- NA
    pur.tolerance <- 0.1
    
    if(match(each.gene, sample.cn.profile$Gene) == 1){
      for(each.pur in c(1:length(pur.vals))){
        if(!is.na(pur.vals[each.pur])) {
          purity.val <- pur.vals[each.pur]
          cn.val <- ploidy.models$Pt[each.pur]
          b.val <- ploidy.models$Bt[each.pur]
          all.models <- append(all.models, list(data.frame("gene"=each.gene, 
                                                           "purity"=as.numeric(purity.val),
                                                           "bval"=as.integer(b.val),
                                                           "cn"=as.integer(cn.val))))
        }
      }
      if(all(is.na(all.models))) next
      names(all.models) <- paste0("model_", c(1:length(all.models)))
    } else {
      all.models <- lapply(names(all.models), function(i.id){
        #i <- all.models[[1]]
        i <- all.models[[i.id]]
        lo.pur <- i$purity - pur.tolerance
        hi.pur <- i$purity + pur.tolerance
        
        pur.vals.idx <- which(pur.vals > lo.pur & pur.vals < hi.pur)
        purity.val <- pur.vals[pur.vals.idx]
        cn.val <- ploidy.models$Pt[pur.vals.idx]
        b.val <- ploidy.models$Bt[pur.vals.idx]
        rbind(all.models[[i.id]], data.frame("gene"=rep(each.gene, length(purity.val)),
                                             "purity"=as.numeric(purity.val),
                                             "bval"=as.integer(b.val),
                                             "cn"=as.integer(cn.val)))

      })
      if(all(is.na(all.models))) next
      names(all.models) <- paste0("model_", c(1:length(all.models)))
    }
  }
  
  # Identify models that fit copy-number restrictions
  bf.model <- lapply(all.models, function(model.x){
    # User specified relative copy-state
    user.cn <- unique(sort(sample.cn.profile$CN))
    mode.model.x <- lapply(user.cn, function(cn.state){
      # Identifies genes with the given relative copy-state
      idx <- which(sample.cn.profile$CN %in% cn.state) 
      gene.idx <- which(model.x$gene %in% sample.cn.profile[idx,]$Gene)
      model.cnstate.x <- model.x[gene.idx,]
      
      # For these genes, identifies the most prevalent ploidy/cn-model that is held constant and selects for these purity/ploidy models
      mode.cn <- getModeCnState(model.cnstate.x)
      mode.cn.idx <- sapply(names(mode.cn), function(nm) which(with(model.cnstate.x, cn == nm)))
      mode.model.x <- tryCatch({ 
        mode.model.x <- model.cnstate.x[mode.cn.idx, ]
        mode.model.x$rank <- cn.state
        return(mode.model.x)
        },error=function(e){as.data.frame(matrix(nrow=0, ncol=5))})
      
      return(mode.model.x)
    })
    
    names(mode.model.x) <- user.cn
    
    # Hack-y sort of a filter for impossible cases
    if(all(sapply(mode.model.x, function(i) nrow(i) > 0))){
      if(length(mode.model.x) > 1){
        min.mode <- unique(mode.model.x[[1]]$cn)
        max.mode <- unique(mode.model.x[[length(user.cn)]]$cn)
        if(min(min.mode) > max(max.mode)) {
          mode.model.x <- list()
        } else {
          rm.idx <- which(mode.model.x[[length(user.cn)]]$cn <= max(min.mode))
          rm.idx <- c(rm.idx, which(mode.model.x[[1]]$cn >= min(max.mode)))
          if(length(rm.idx) > 0) mode.model.x[[length(user.cn)]] <- mode.model.x[[length(user.cn)]][-rm.idx,]
        }
      }
    }
    
    do.call("rbind", mode.model.x)
  })
  
  # Remove models that have 1 or more genes that don't fit the model
  if(length(bf.model) > 1){
    incomplete.models <- which(sapply(bf.model, function(x) !all(sample.cn.profile$Gene %in% x$gene)))
    
    # return(list("complete" = bf.model[-incomplete.models],
    #             "incomplete" = bf.))
    best.fit.models[['complete']][[net.id]] <- if(length(incomplete.models) > 0) bf.model[-incomplete.models] else bf.model
    best.fit.models[['incomplete']][[net.id]] <- bf.model
  }


}



