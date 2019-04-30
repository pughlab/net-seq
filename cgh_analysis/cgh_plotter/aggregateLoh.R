library(taRifx)
library(scales)
library(ltm)
#### FUNCTIONS ####
mutColours <- function(){
 list("INS"="#82CABC",
      "DEL"="#F56A68",
      "NON"="#157F7D",
      "MIS"="#72A3C5",
      "SPLICE"="#F1A85B",
      "START"="#F1A85B",
      "NA"="white") 
}

lohColours <- function(){
  c("0"=alpha('white', 0),
    "1"='dodgerblue',
    "2"='thistle',
    "3"='red',
    "4"='#489C47')
}

missegregateChr <- function(){
  ## Transcribed and estimated using Inkscape Scaling
  ## directly from Figure 1.I of the Worrell Paper
  ## representing % of Cells for single-cell nocodazole treatment
  chrcols <- c(4.75, 2.8, 2.8, 0.8, 0.8, 2.1,
               0, 2.1, 0.8, 1.5, 2.1, 1.5,
               2.1, 0, 0, 0, 0.8, 2.8, 0.8,
               2.1, 0.8, 0.8, 0.8)
  names(chrcols) <- c(chrs, "chrX")
  
  chrcols <- round(chrcols/max(chrcols), 3)
  chrcols
}

lohChr <- function(){
  loh=paste0("chr", c(1, 2, 3, 6, 8, 10, 11, 16, 18, 21, 22))
  ret=paste0("chr", c(4, 5, 7, 9, 12, 13, 14, 15, 17, 19, 20))

  chrcols <- rep(NA, length(chrs))
  names(chrcols) <- chrs
  chrcols[match(loh, chrs)] <- 'black'
  chrcols[match(ret, chrs)] <- 'white'
  chrcols
}

toFactor <- function(x){
  as.numeric(as.character(x))
}

plotBlank <- function(chr, xlim=c(0,1), ylim=NULL){
  if(is.null(ylim)) ylim <- (-1 * c(ucsc.spl[[chr]]$size, 1))
  plot(0, type='n', 
       xlim=xlim,  ylim=ylim, 
       axes=FALSE, lty=0, ylab='', xlab='')
}

plotChr <- function(chr){
  plotBlank(chr)
  loh.tmp <- remove.factors(seg.spl[[chr]][which(seg.spl[[chr]]$LOH == 1),])
  if(is.null(loh.tmp)) loh.tmp <- matrix(nrow = 0, ncol=1)
  if(nrow(loh.tmp) > 0){
    if(any(loh.tmp$CNt > 3) && each.disease != 'ccl') loh.tmp[which(loh.tmp$CNt > 3),]$CNt <- 3
    rect(xleft = 0, ybottom = (-1 * as.integer(loh.tmp$Start.bp)), 
         xright = 1, ytop = (-1 * as.integer(loh.tmp$End.bp)), 
         border=NA, col=lohColours()[as.character(loh.tmp$CNt)])
  }
}

plotChrAxis <- function(chr){
  plot(0, type='n', 
       xlim=c(0,10),  ylim=c(0, 10.5), 
       axes=FALSE, lty=0, ylab='', xlab='')
  text(x=10, y = 5, labels = gsub("^chr", "", chr, ignore.case = TRUE), 
       pos = 2, cex=chr.cex)
  if(chr=='chr1') segments(x0 = 9,y0 = 10,x1 = 10,y1 = 10)
  segments(x0 = 9,y0 = 0,x1 = 10,y1 = 0)
  
  if(each.disease == 'pnets'){
    if(chr=='chr1') text(x=c(2,4), y=c(9.5, 9.5), labels=c("L", "M"), cex=0.5)
    points(x=4, y=4, bg=alpha("black", missegregateChr()[chr]), pch=21)
    points(x=2, y=4, bg=lohChr()[chr], pch=21)
  }
}

getMargins <- function(disease.idx){
  sample.id <- names(disease.list[[each.disease]])[disease.idx]
  if(length(sample.id) == 0) sample.id <- 'x'
  if(!is.na(match(sample.id, names(grp)))){
    p <- switch(grp[[sample.id]],
                "a"=c(0, sp1, 0, sp),
                "b"=c(0, sp, 0, sp1))
  } else {
    p <- c(0,sp,0,sp)
  }
  p
}

addRect <- function(inds, cols, ylim=c(0,2)){
  plotBlank('chr1', xlim=c(0,1), ylim=ylim)
  sapply(seq_along(inds), function(i){
    idx <- inds[[i]]
    rect(xleft = 0, ybottom = idx[1], 
         xright = 1, ytop = idx[2], 
         col=cols[i], border=FALSE)
  })
}

fixSeg <- function(seg){
  seg$modal_A1 <- rmFactor(seg$modal_A1)
  seg$Chromosome <- as.character(seg$Chromosome)
  seg$LOH <- rmFactor(seg$LOH)
  
  #modal_A2
  if(all(is.na(seg$modal_A2))){
    seg$modal_A2 <- 1
  }
  seg$modal_A2 <- rmFactor(seg$modal_A2)
  
  # Chromosome
  if(!any(grepl("chr", seg$Chromosome))){
    seg$Chromosome <- paste0("chr", seg$Chromosome)
  }
  seg <- seg[order(factor(seg$Chromosome, levels=chrs)),]
  
  # length
  seg$End.bp <- rmFactor(seg$End.bp) 
  seg$Start.bp <- rmFactor(seg$Start.bp)
  seg$length <- with(seg, End.bp - Start.bp)
  
  # CNt
  if(!any(grepl("CNt", colnames(seg)))){
    seg$modal_A1 <- rmFactor(seg$modal_A1) 
    seg$modal_A2 <- rmFactor(seg$modal_A2)
    seg$CNt <- with(seg, modal_A2 + modal_A1)
  }
  if(each.disease == 'ccl') seg$CNt <- '4'

  seg
}

rmFactor <- function(x){
  as.numeric(as.character(x))
}


#### VARIABLES ####
mut.table <- read.csv("~/git/net-seq/cgh_analysis/data/pnet_mutations.csv",
                      header=TRUE, stringsAsFactors = FALSE, check.names = FALSE)
load("~/git/net-seq/cgh_analysis/data/chromInfo.Rdata")
load("~/git/net-seq/cgh_analysis/data/data.aggregateLoh.pnets.Rdata")

ucsc.spl <- split(ucsc.chrom, f=ucsc.chrom$chrom)
chrs <- paste0("chr", c(1:22))
sp <- 0.005
sp1 <- 0.5
chr.cex <- 0.7

bisir <- biserial.cor(missegregateChr()[-23], as.integer(factor(lohChr())))

disease.names <- c("pnets", "otb-pnet", "ccl")
for(each.disease in disease.names){
  pdf(paste0("~/Desktop/netseq/loh-cn/", each.disease, "_loh.pdf"), 
      height = 7, width = 2 + (0.3 * length(names(disease.list[[each.disease]]))))


  if(each.disease == 'pnets'){
    grp <- list("NET-003-T_1"="a",
                "NET-003-T_2"="b",
                "NET-009-T_1"="a",
                "NET-009-T_2"="b")
    
    sample.ord <- c('008', '001', '003', '009')
    sample.ord <- unlist(sapply(sample.ord, function(i) {
      grep(i, names(disease.list[[each.disease]]))
    }))
    disease.list[[each.disease]] <- disease.list[[each.disease]][sample.ord]
  } else if(each.disease == 'otb-pnet'){
    keep.idx <- names(disease.list[[each.disease]]) %in% colnames(mut.table)
    disease.list[[each.disease]] <- disease.list[[each.disease]][keep.idx]
    
    sample.ord <- order(factor(names(disease.list[[each.disease]]), 
                               levels=colnames(mut.table)))
    disease.list[[each.disease]] <- disease.list[[each.disease]][sample.ord]
  } else if(each.disease =='ccl'){
    names(disease.list[[each.disease]]) <- gsub("(?<!-)1$", "-1", 
                                                names(disease.list[[each.disease]]), perl=TRUE)
    
  }
  
  #### VISUALIZATION ####
  {
    close.screen(all.screens=TRUE)
    main.spl <- split.screen(matrix(c(0,1,0.85,1,
                                      0,1,0.75,0.85,
                                      0,1,0.2,0.75,
                                      0,1,0,0.2), byrow = TRUE, ncol=4))
    mut.screen <- main.spl[1]
    frc.screen <- main.spl[2]; gen.loh <- list()
    loh.screen <- main.spl[3]
    sid.screen <- main.spl[4]
    
    #### Mutation ####
    ## Plot the Mutation track
  
    screen(mut.screen)
    dis.scrn <- split.screen(c(1, (length(disease.list[[each.disease]]) + 1)))
    for(disease.idx in dis.scrn){
      screen(disease.idx);
      
      disease.idx <- (disease.idx - min(dis.scrn))
      p <- getMargins(disease.idx)
      p[c(1, 3)] <- c(0.5, 1) # add bottom padding
      par(mar=p)
      
      if(disease.idx == 0){
        ## Plot the Diagnosis/Second biopsy labels
        plotBlank('chr1', xlim=c(0,10), ylim=c(0,4))
        text(x=rep(10, 4), y = c(0.5, 1.5, 2.5, 3.5), 
             cex=0.5, pos = 2, font=3,
             labels = c("TP53", "ATRX", "DAXX", "MEN1"))
      } else {
        ## Plot group inclusion
        sid <- names(disease.list[[each.disease]])[disease.idx]
        muts <- mut.table[,match(sid, colnames(mut.table))]
        muts[is.na(muts)] <- "NA"
        idx <- list(c(0,1), c(1,2), c(2,3), c(3,4))
        
        addRect(idx, cols = rev(unlist(mutColours()[muts])), ylim=c(0,4))
        axis(side = 1, at=c(0,1), lwd.ticks = 0, labels=c("", ""))
      }
    }
    
    
    #### Sample ID ####
    ## Plot the Sample ID track
    screen(sid.screen)
    dis.scrn <- split.screen(c(1, (length(disease.list[[each.disease]]) + 1)))
    for(disease.idx in dis.scrn){
      screen(disease.idx);
      
      disease.idx <- (disease.idx - min(dis.scrn))
      p <- getMargins(disease.idx)
      p[c(1, 3)] <- c(4.1, 0.5) # add bottom padding
      par(mar=p)
      
      if(disease.idx == 0 & each.disease == 'pnets'){
        ## Plot the Diagnosis/Second biopsy labels
        plotBlank('chr1', xlim=c(0,10), ylim=c(0,2))
        text(rep(10, 2), y = c(0.5, 1.5), cex=0.5, pos = 2,
             labels = c("2nd biopsy", "Diagnosis"))
      } else if(disease.idx != 0) {
        ## Plot group inclusion
        sid <- names(disease.list[[each.disease]])[disease.idx]
        sid <- strsplit(sid, split="-")[[1]]
        idx <- switch(sid[3],
                      "T_1" = c(1, 2),
                      "T_2" = c(0, 1))
        if(each.disease == 'pnets') {
          addRect(list(idx), "black", ylim=c(0,2))
          axis(side = 3, at=c(0,1), lwd.ticks = 0, labels=c("", ""))
        } 
        if(each.disease == 'otb-pnet') sid[2] <- gsub("^0", "1", sid[2])
        axis(side = 1, at = 0.5, labels = paste(sid[1:2], collapse="-"), 
             las=2, cex.axis=0.7)
      }
    }
    
    #### LOH ####
    ## Plot the LOH Track
    screen(loh.screen)
    dis.scrn <- split.screen(c(1, (length(disease.list[[each.disease]]) + 1)))
    for(disease.idx in dis.scrn){
      screen(disease.idx);
      
      disease.idx <- (disease.idx - min(dis.scrn))
      p <- getMargins(disease.idx)
      
      if(disease.idx == 0){
        ## Plot the chromosome axes
        chr.scrn <- split.screen(c(length(chrs), 1))
        lapply(chrs, function(chr){
          screen(chr.scrn[match(chr, chrs)]); par(mar=c(0, 0, 0, 0))
          plotChrAxis(chr)
        })
      } else {
        ## Plot individual sample LOH per chromosome
        seg <- disease.list[[each.disease]][[disease.idx]]
        seg <- fixSeg(seg)
        
        gen.loh[[disease.idx]] <- sum(seg[which(seg$LOH==1),]$length) / total.genome.size
        seg.spl <- split(seg, f=seg$Chromosome)
        
        chr.scrn <- split.screen(c(length(chrs), 1))
        lapply(chrs, function(chr){
          screen(chr.scrn[match(chr, chrs)]); par(mar=p)
          plotChr(chr)
        })
      }
    }
    
    #### Fraction LOH ####
    ## Plot the Genomic LOH Track
    screen(frc.screen)
    dis.scrn <- split.screen(c(1, (length(disease.list[[each.disease]]) + 1)))
    for(disease.idx in dis.scrn){
      screen(disease.idx);
      
      disease.idx <- (disease.idx - min(dis.scrn))
      p <- getMargins(disease.idx)
      p[c(1,3)] <- c(0.5,0.1)
      par(mar=p)
      
      if(disease.idx == 0){
        ## Plot the Axis labels for the barplot
        plotBlank('chr1', xlim=c(0,10), ylim=c(0,11))
        gen.frac.scale <- seq(0, 10, by=2)
        text(x=rep(9, length(gen.frac.scale)), y=gen.frac.scale, 
             labels = gen.frac.scale*(100/max(gen.frac.scale)), 
             pos = 2, cex=0.5)
        segments(x0 = rep(9, length(gen.frac.scale)), y0 = gen.frac.scale ,
                 x1 = rep(10, length(gen.frac.scale)), y1 = gen.frac.scale)
      } else {
        ## Plot group inclusion
        plotBlank(chr, xlim=c(0,1), ylim=c(0,1.1))
        rect(xleft = 0, ybottom = 0,xright = 1,ytop = gen.loh[[disease.idx]],
             col="gray40", border=FALSE)
      }
    }
    close.screen(all.screens=TRUE)
  }
  dev.off()
}
