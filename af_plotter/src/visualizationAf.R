# Function: plotAfExomeRna
# Purpose:  Identical to plotAf, but will plot RNA and Exome data on the same graph for comparison
# Input:  name: Name of the sample/output file
#         exome.list: The returned list from generateChrList for exome data
#         rna.list: The returned list from generateChrList for rnaseq data
#       max.snp.density = max snps per kb (default = 5)
#       min.snp.density = min snps per kb (default = 5)
# Returns:  Output plots to the output file
plotAfExomeRna <- function(name="unlabelled", exome.list=list(), rna.list=list(),
                           ebin.list=NULL, rbin.list=NULL,
                           min.depth=15, filt.region=FALSE, 
                           plot.bins=FALSE, plot.points=TRUE,
                           max.snp.density=5, min.snp.density=0.05, bin.size=50000,
                           plot.type='png'){
  ####################
  #### OUTPUT: Plot all Chr together
  require(scales)
  if(plot.type=='png'){
    png(file=paste(name, ".min", paste=min.depth, ".allChr.png", sep=""), width=2000, height=1750, res=200)
  } else if(plot.type=='pdf'){
    pdf(file=paste(name, ".min", paste=min.depth, ".allChr.pdf", sep=""))
  }
  
  layout(matrix(c(seq(1:(length(exome.list) + 2)), 
                  (seq(1:(length(exome.list) + 2)) + (length(exome.list) + 2))), nrow=2, byrow=T))
  rna.exome.list <- list(exome.list, rna.list)
  if(rm.x){
    chr.num <- paste("chr", c(1:22), sep="")
  } else {
    chr.num <- paste("chr", c(1:22, "X"), sep="")
  }
  
  for(z.af.list in rna.exome.list){
    par(mar=c(3,2.1,3,0))
    plot(c(0), xaxt='n', ylim=c(0,1))
    chr.col.cnt <- 1
    for(each.chr in chr.num){
      each.chr.df <- z.af.list[[each.chr]]
      par(mar=c(3,0,3,0))
      
      each.chr.df <- fixFactors(each.chr.df, c("depth", "pos"), 'integer')
      each.chr.df <- fixFactors(each.chr.df, c("a", "b"), 'numeric')
      conf.chr.df <- each.chr.df[which(each.chr.df$depth > min.depth),]
      
      plot(x=conf.chr.df$pos,
           y=conf.chr.df$a,
           ylim=c(0,1),
           col=alpha("black", alpha=0.5),
           pch=16,
           yaxt='n', xaxt='n')
      points(x=conf.chr.df$pos,
             y=conf.chr.df$b,
             pch=16,
             col=alpha("grey", alpha=0.5))
      
      rect(par("usr")[1],par("usr")[3],par("usr")[2],0,col = alpha("grey", 0.6))
      if(chr.col.cnt %% 2 == 1){
        rect(par("usr")[1],par("usr")[3],par("usr")[2],par("usr")[4],col = alpha("grey", 0.2))    
        
        text(x=(((par("usr")[2] - par("usr")[1]) / 2) + par("usr")[1]), 
             y=((par("usr")[3] - 0) / 2),
             labels=paste("Chr", chr.col.cnt, sep=""),
             cex=0.8)
      }
      
      # Add lines indicating copy ratio
      abline(h=c(0.2, 0.33, 0.5, 0.66, 0.8), col="grey")
      chr.col.cnt <- chr.col.cnt + 1
    }
    par(mar=c(3,0,3,0))
    plot.new()
  }
  dev.off()
  
  ####################
  #### OUTPUT: Plot Individual Chr
  dna.seg <- list() # List contain all AF segments called by CBS for DNA
  rna.seg <- list() # List contain all AF segments called by CBS for RNA
  chr.col.cnt <- 1
  for(each.chr in chr.num){
    if(plot.type=='png'){
      png(file=paste(name, ".min", paste=min.depth, ".", each.chr, ".png", sep=""), 
          width=2000, height=1750, res=200)
    } else if(plot.type=='pdf'){
      pdf(file=paste(name, ".min", paste=min.depth, ".", each.chr, ".pdf", sep=""))
    }
    
    layout(matrix(c(0,0,4,4,4,4,9,9,9,9,0,0,
                    0,0,4,4,4,4,9,9,9,9,0,0,
                    0,0,3,3,3,3,8,8,8,8,0,0,
                    0,2,1,1,1,1,6,6,6,6,7,0,
                    0,2,1,1,1,1,6,6,6,6,7,0,
                    0,2,1,1,1,1,6,6,6,6,7,0,
                    0,2,1,1,1,1,6,6,6,6,7,0,
                    0,2,1,1,1,1,6,6,6,6,7,0,
                    0,2,1,1,1,1,6,6,6,6,7,0,
                    0,2,1,1,1,1,6,6,6,6,7,0,
                    0,2,1,1,1,1,6,6,6,6,7,0,
                    0,2,1,1,1,1,6,6,6,6,7,0,
                    0,0,5,5,5,5,10,10,10,10,0,0,
                    0,0,5,5,5,5,10,10,10,10,0,0,
                    0,0,5,5,5,5,10,10,10,10,0,0,
                    0,0,14,14,14,14,19,19,19,19,0,0,
                    0,0,14,14,14,14,19,19,19,19,0,0,
                    0,0,13,13,13,13,18,18,18,18,0,0,
                    0,12,11,11,11,11,16,16,16,16,17,0,
                    0,12,11,11,11,11,16,16,16,16,17,0,
                    0,12,11,11,11,11,16,16,16,16,17,0,
                    0,12,11,11,11,11,16,16,16,16,17,0,
                    0,12,11,11,11,11,16,16,16,16,17,0,
                    0,12,11,11,11,11,16,16,16,16,17,0,
                    0,12,11,11,11,11,16,16,16,16,17,0,
                    0,12,11,11,11,11,16,16,16,16,17,0,
                    0,12,11,11,11,11,16,16,16,16,17,0,
                    0,0,15,15,15,15,20,20,20,20,0,0,
                    0,0,15,15,15,15,20,20,20,20,0,0,
                    0,0,15,15,15,15,20,20,20,20,0,0), nrow=30, byrow=T))
    
    for(each.data in c("dna", "rna")){
      # Subset the dna list for each Chromosome and Arm
      each.chr.df <- list()
      if(each.data %in% 'dna'){
        each.chr.df <- exome.list[[each.chr]]
        each.chr.bin.df <- ebin.list[[each.chr]]
        ebin.list[[each.chr]]$arm <- NA
      } else if(each.data %in% 'rna'){
        each.chr.df <- rna.list[[each.chr]]
        each.chr.bin.df <- rbin.list[[each.chr]]
        rbin.list[[each.chr]]$arm <- NA
      }
      
      #Change factors to integer/numeric
      each.chr.df <- fixFactors(each.chr.df, c("depth", "pos"), 'integer')
      each.chr.df <- fixFactors(each.chr.df, c("a", "b", "raw_a", "raw_b"), 'numeric')
      
      # Removes variant SNPs below a certain threshold
      conf.chr.df <- mergeAf(each.chr.df, NULL, min.depth)
      
      # Filters out regions where there are a LOT of SNPs in a 5k bin
      if(filt.region == TRUE & plot.bins == TRUE){
        #SNP Density: http://bmcgenomics.biomedcentral.com/articles/10.1186/1471-2164-8-146    Fig. 2
        min.bin.n <- min.snp.density * (bin.size / 1000)
        max.bin.n <- max.snp.density * (bin.size / 1000)
        each.chr.bin.df <- each.chr.bin.df[which(each.chr.bin.df$n > min.bin.n &
                                                   each.chr.bin.df$n < max.bin.n),]
        #conf.chr.df <- filtRegion(conf.chr.df, bin=5000)
      }
#       max.depth <- quantile(as.integer(each.chr.df$depth), 0.90)
#       #conf.chr.df <- conf.chr.df[which(conf.chr.df$depth > min.depth),]
#       conf.chr.df <- conf.chr.df[which(conf.chr.df$depth > max.depth),]
      
      chr.arms.df <- chrom.df[which(chrom.df$chr == each.chr),]
      
      seg.pos.a <- c()
      seg.pos.b <- c()
      snp.pos.tot <- c()
      #Plot Arm Level Plots
      for(each.arm in c("p", "q")){
        spos <- chr.arms.df[which(chr.arms.df$arm == each.arm), 'start_pos']
        epos <- chr.arms.df[which(chr.arms.df$arm == each.arm), 'end_pos']
        
        #retrieve snps on the chromosomal arm
        
        if(plot.bins){
          bin.a <- each.chr.bin.df[which(each.chr.bin.df$start.pos >= spos & 
                                         each.chr.bin.df$end.pos <= epos), 'mean.a']
          bin.b <- each.chr.bin.df[which(each.chr.bin.df$start.pos >= spos & 
                                         each.chr.bin.df$end.pos <= epos), 'mean.b']
          bin.pos <- each.chr.bin.df[which(each.chr.bin.df$start.pos >= spos & 
                                             each.chr.bin.df$end.pos <= epos), 'start.pos']
        }
        snp.a <- conf.chr.df[which(conf.chr.df$pos >= spos & conf.chr.df$pos <= epos), 'a']
        snp.b <- conf.chr.df[which(conf.chr.df$pos >= spos & conf.chr.df$pos <= epos), 'b']
        snp.pos <- conf.chr.df[which(conf.chr.df$pos >= spos & conf.chr.df$pos <= epos), 'pos']
        
        if(length(snp.pos) > 1){
          # Retrieve the pathogenic SNPs from clinvar annotation
          risk.match <- c("Pathogenic", "risk")
          risk.rows <- grep(paste(risk.match,collapse="|"), conf.chr.df$clinvar_state)
          gene.pos.start <- 0; gene.pos.end <- 0; gene.pos.name <- ""; gene.name.y <- 0
          if(length(risk.rows) > 0){
            risk.rows <- risk.rows[risk.rows <= max(which(conf.chr.df$pos <= epos)) &
                                     risk.rows >= min(which(conf.chr.df$pos >= spos))]
            risk.rows <- risk.rows[!is.na(conf.chr.df[risk.rows,'clinvar_gene_start'])]
            if(length(risk.rows) > 0){
              gene.pos.start <- conf.chr.df[risk.rows,'clinvar_gene_start']
              gene.pos.end <- conf.chr.df[risk.rows, 'clinvar_gene_end']
              gene.pos.name <- conf.chr.df[risk.rows, 'clinvar_gene']
              gene.name.y <- as.vector(apply(matrix(c(0.3, 0.6), nr=2), 2, rep, length(risk.rows)/2)) 
              if(length(risk.rows) %% 2 == 1) gene.name.y <- c(gene.name.y, 0.3) else gene.name.y
            } 
          }        
          
          ####################
          #### Plotting the AF for each chr
          if(each.arm == 'p'){
            par(mar=c(0, 0, 0, 0))
          } else {
            par(mar=c(0, 0, 0, 0))
          }        
          # plot SNP points
          plot(0, type='n', ylim=c(0,1), xlim=c(spos, epos),
               xaxt='n', yaxt='n')
          if(plot.bins){
            rect(xleft = each.chr.bin.df$start.pos, ybottom = each.chr.bin.df$low.sd.a, 
                 xright = each.chr.bin.df$end.pos, ytop = each.chr.bin.df$high.sd.a, 
                 col=alpha("cyan", 0.3), lty=0)
            rect(xleft = each.chr.bin.df$start.pos, ybottom = (bin.a - 0.01), 
                 xright = each.chr.bin.df$end.pos, ytop = (bin.a + 0.01), 
                 col="black", lty=0)
            rect(xleft = each.chr.bin.df$start.pos, ybottom = each.chr.bin.df$low.sd.b, 
                 xright = each.chr.bin.df$end.pos, ytop = each.chr.bin.df$high.sd.b, 
                 col=alpha("cyan", 0.3), lty=0)
            rect(xleft = each.chr.bin.df$start.pos, ybottom = (bin.b - 0.01), 
                 xright = each.chr.bin.df$end.pos, ytop = (bin.b + 0.01),
                 col="black", lty=0)
            
            seg.list.a <- c(seg.pos.a, plotCbsSegments(each.chr, bin.a, bin.pos, each.data, 'a'))
            seg.list.b <- c(seg.pos.b, plotCbsSegments(each.chr, bin.b, bin.pos, each.data, 'b'))
            snp.pos.tot <- c(snp.pos.tot, bin.pos)
          }
          if(plot.points){
            points(x=snp.pos,
                   y=snp.a,
                   col=alpha("grey", 0.5))
            points(x=snp.pos,
                   y=snp.b,
                   col=alpha("grey", 0.5))
            seg.list.a <- c(seg.pos.a, plotCbsSegments(each.chr, snp.a, snp.pos, each.data, 'a'))
            seg.list.b <- c(seg.pos.b, plotCbsSegments(each.chr, snp.b, snp.pos, each.data, 'b'))
            snp.pos.tot <- c(snp.pos.tot, snp.pos)
          }
          seg.pos.a <- c(seg.pos.a, seg.list.a$seg.pos)
          seg.pos.b <- c(seg.pos.b, seg.list.b$seg.pos)
          
          # plots Pathogenic SNP points
          points(x=conf.chr.df[grep("Pathogenic", conf.chr.df$clinvar_state), 'pos'],
                 y=conf.chr.df[grep("Pathogenic", conf.chr.df$clinvar_state), 'b'],
                 col="blue")
          points(x=conf.chr.df[grep("Pathogenic", conf.chr.df$clinvar_state), 'pos'],
                 y=conf.chr.df[grep("Pathogenic", conf.chr.df$clinvar_state), 'a'],
                 col="blue")
          # plots SNV points
          points(x=conf.chr.df[which(!is.na(conf.chr.df$Variant_Classification)), 'pos'],
                 y=conf.chr.df[which(!is.na(conf.chr.df$Variant_Classification)), 'a'],
                 col="red")
          points(x=conf.chr.df[which(!is.na(conf.chr.df$Variant_Classification)), 'pos'],
                 y=conf.chr.df[which(!is.na(conf.chr.df$Variant_Classification)), 'b'],
                 col="red")
          abline(h=c(0.2, 0.33, 0.5, 0.66, 0.8), col="grey")

          
          
          # Append the DNA segments from CBS to the end of the return list
          #seg.list.a$segment.smoothed$output <- cbind(seg.list.a$segment.smoothed$output, seg.list.b$segment.smoothed$output$seg.mean)
          if(!is.null(seg.list.a$segment.smoothed$output)){
            seg.list.a$segment.smoothed$output <- cbind(seg.list.a$segment.smoothed$output, (1-seg.list.a$segment.smoothed$output$seg.mean))
            colnames(seg.list.a$segment.smoothed$output) <- c("ID", "chrom", "loc.start", "loc.end", "num.mark", "seg.mean.a", "seg.mean.b")
            if(each.data %in% 'dna'){
              dna.seg[[each.chr]] <- rbind(dna.seg[[each.chr]], seg.list.a$segment.smoothed$output)
            } else if(each.data %in% 'rna'){
              rna.seg[[each.chr]] <- rbind(rna.seg[[each.chr]], seg.list.a$segment.smoothed$output)
            }
          }
          
          
          ####################
          #### Plot density of the AF
          if(plot.bins){
            dens.a <- density(bin.a, na.rm=TRUE)
            dens.b <- density(bin.b, na.rm=TRUE)
          }
          if(plot.points){
            dens.a <- density(snp.a, na.rm=TRUE)
            dens.b <- density(snp.b, na.rm=TRUE)
          }
          
          if(each.arm == 'p'){
            par(mar=c(0, 0, 0, 0))
          } else {
            par(mar=c(0, 0, 0, 0))
          }
          plot(if(each.arm == 'q') dens.a$y else -dens.a$y, 
               dens.a$x, 
               type="l", col="black",
               ylab=if(each.arm == 'q') 'Density' else '',
               ylim=c(0,1),
               yaxt='n', xaxt='n')
          lines(if(each.arm == 'q') dens.b$y else -dens.b$y, 
                dens.b$x, 
                col="black")
          # Lines indicating where the 1:1, 2:1, 1:0 ratios are
          abline(h=c(0.2, 0.33, 0.5, 0.66, 0.8), col="grey")
          if(each.arm == 'q'){
            axis(4, at=c(0.20, 0.33, 0.50, 0.66, 0.80),
                 labels=c(0.20, 0.33, 0.50, 0.66, 0.80), las=2)
            
            # Density label on x-axis
            axis(1, at=3, labels="Density", tick=F, )
          } else if(each.arm == 'p'){
            axis(2, at=c(0, 0.2, 0.4, 0.6, 0.8, 1.0),
                 labels=c(0, 0.2, 0.4, 0.6, 0.8, 1.0), las=2)
          }
          
          ####################
          ####Track indicating location of SNV Genes
          par(mar=c(0, 0, 0, 0))
          plot(0, 0, type='n',
               xlim=c(spos, epos),
               ylim=c(0,1),
               yaxt='n', xaxt='n',
               ylab="SNV", las=2)
          if(dim(conf.chr.df[which(!is.na(conf.chr.df$Variant_Classification)),])[1] > 0){
            rect(xleft=(conf.chr.df[which(!is.na(conf.chr.df$Variant_Classification)), 'pos'] - 20000), 
                 ybottom=0, 
                 xright=(conf.chr.df[which(!is.na(conf.chr.df$Variant_Classification)), 'pos'] + 20000), 
                 ytop=1,
                 col=alpha("red", 0.5), border=alpha("red", 0.5))
            snv.name.y <- as.vector(apply(matrix(c(0.3, 0.6), nr=2), 2, rep, length(conf.chr.df[which(!is.na(conf.chr.df$Variant_Classification)), 'pos'])/2))
            if(length(conf.chr.df[which(!is.na(conf.chr.df$Variant_Classification)), 'pos']) %% 2 == 1){
              snv.name.y <- c(snv.name.y, 0.3)
            }
            text(x = conf.chr.df[which(!is.na(conf.chr.df$Variant_Classification)), 'pos'],
                 y = snv.name.y,
                 labels=conf.chr.df[which(!is.na(conf.chr.df$Variant_Classification)), 'Hugo_Symbol'],
                 col="black",
                 cex=0.5)
          }
          
          
          ####################
          ####Track indicating location of Pathogenic Genes
          if(each.arm == 'p'){
            par(mar=c(0, 0, 2.1, 0))
          } else {
            par(mar=c(0, 0, 2.1, 0))
          }  
          plot(0, 0, type='n',
               xlim=c(spos, epos),
               ylim=c(0,1),
               main=paste(if(each.arm=='p') paste(toupper(each.data), ': ', sep="") else '',
                          "Chr", chr.col.cnt, each.arm, sep=""),
               yaxt='n', xaxt='n',
               ylab="SNP", las=2)
          if(!any(gene.pos.start == 0)){
            rect(xleft=gene.pos.start, 
                 ybottom=0, 
                 xright=gene.pos.end, 
                 ytop=1,
                 col=alpha("blue", 0.5), border=alpha("blue", 0.5))
            text(x = (gene.pos.start + ((gene.pos.end - gene.pos.start) * 0.1)),
                 y = gene.name.y,
                 labels=gene.pos.name,
                 col="black",
                 cex=0.5)
          }
          
          ####################
          #### Barplot for depth of each count
          if(each.arm == 'p'){
            par(mar=c(3.1, 0, 0, 0))
          } else {
            par(mar=c(3.1, 0, 0, 0))
          }
          x    <- snp.pos                       # x-coordinates
          y    <- conf.chr.df[which(conf.chr.df$pos >= spos & conf.chr.df$pos <= epos), 'depth']                                 # y-coordinates
          barW <- 1 / max(snp.pos)                                   # width of bars
          plot(x, y, type='n',
               ylab=if(each.arm=='q') '' else 'Depth',
               xlab="Genomic Position",
               xlim=c(spos, epos),
               ylim=c(0,400),
               yaxt='n', xaxt='n')
          rect(xleft=x-barW, 
               ybottom=0, 
               xright=x+barW, 
               ytop=conf.chr.df$depth,
               col=gray(0.5))
          lines(x, y, col=alpha("red", 0.4))
  
        } else {
          #######Create Empty Plots
          if(each.arm == 'p'){
            par(mar=c(0, 0, 0, 0))
          } else {
            par(mar=c(0, 0, 0, 0))
          }   
          plot(1, ylim=c(0,1),
               xlim=c(spos, epos),
               ylab=if(each.arm == "q") "" else 'Allelic Fraction',
               xaxt='n',
               yaxt=if(each.arm == "q") 'n' else 't')
          abline(h=c(0.2, 0.33, 0.5, 0.66, 0.8), col="grey")
          
          # Plot density of the AF
          if(each.arm == 'p'){
            par(mar=c(0, 0, 2.1, 0))
          } else {
            par(mar=c(0, 0, 2.1, 0))
          }
          plot(1,ylab=if(each.arm == 'q') 'Density' else '',
               ylim=c(0,1),
               yaxt='n', xaxt='n')
          # Lines indicating where the 1:1, 2:1, 1:0 ratios are
          abline(h=c(0.2, 0.33, 0.5, 0.66, 0.8), col="grey")
          if(each.arm == 'q'){
            axis(4, at=c(0.20, 0.33, 0.50, 0.66, 0.80),
                 labels=c(0.20, 0.33, 0.50, 0.66, 0.80))
            # Density label on x-axis
            if(each.arm == 'q') axis(1, at=1, labels="Density")
          }
          
          # Track indiciating blank SNVs
          par(mar=c(0, 0, 0, 0)) 
          plot(0, 0, type='n',
               xlim=c(spos, epos),
               ylim=c(0,1),
               yaxt='n', xaxt='n')
          
          #Track indicating location of Pathogenic Genes
          if(each.arm == 'p'){
            par(mar=c(0, 0, 2.1, 0))
          } else {
            par(mar=c(0, 0, 2.1, 0))
          }  
          plot(0, 0, type='n',
               xlim=c(spos, epos),
               ylim=c(0,1),
               main=paste(if(each.arm=='p') paste(toupper(each.data), ': ', sep="") else '',
                          "Chr", chr.col.cnt, each.arm, sep=""),
               yaxt='n', xaxt='n')
          
          # Barplot for depth of each count
          if(each.arm == 'p'){
            par(mar=c(3.1, 0, 0, 0))
          } else {
            par(mar=c(3.1, 0, 0, 0))
          }
          plot(1, type='n',
               ylab=if(each.arm=='q') '' else 'Depth',
               xlab="Genomic Position",
               xlim=c(spos, epos),
               ylim=c(0,400),
               yaxt=if(each.arm=='q') 'n' else 't')
        }
        
        # Add Arm information to the dataframe
        if(each.data %in% 'dna' & plot.bins){
          ebin.list[[each.chr]][which(ebin.list[[each.chr]]$start.pos %in% 
                                        bin.pos),'arm'] <- rep(each.arm, length(bin.pos))
        } else if(each.data %in% 'rna' & plot.bins){
          rbin.list[[each.chr]][which(rbin.list[[each.chr]]$start.pos %in%  
                                        each.chr.bin.df$start.pos),'arm'] <- rep(each.arm, dim(each.chr.bin.df)[1])
        }
        
      }
      
      #Add the segment positions called from the CBS algorithm to the list for exome or rna for each chr
      if(each.data %in% 'dna'){
        #exome.list[[each.chr]][which(as.integer(as.character(exome.list[[each.chr]]$depth)) > min.depth),'seg.a'] <- seg.pos.a
        #exome.list[[each.chr]][which(as.integer(as.character(exome.list[[each.chr]]$depth)) > min.depth),'seg.b'] <- seg.pos.b
        if(plot.points){ 
          exome.list[[each.chr]][which(exome.list[[each.chr]]$pos %in% snp.pos.tot),'seg.a'] <- seg.pos.a
          exome.list[[each.chr]][which(exome.list[[each.chr]]$pos %in% snp.pos.tot),'seg.b'] <- seg.pos.b
        }
        if(plot.bins){
          ebin.list[[each.chr]][which(ebin.list[[each.chr]]$start.pos %in% snp.pos.tot), 'seg.a'] <- seg.pos.a
          ebin.list[[each.chr]][which(ebin.list[[each.chr]]$start.pos %in% snp.pos.tot), 'seg.b'] <- seg.pos.b
        }        
      } else if(each.data %in% 'rna'){
        #rna.list[[each.chr]][which(as.integer(as.character(rna.list[[each.chr]]$depth)) > min.depth),'seg.a'] <- seg.pos.a
        #rna.list[[each.chr]][which(as.integer(as.character(rna.list[[each.chr]]$depth)) > min.depth),'seg.b'] <- seg.pos.b
        if(plot.points){ 
          rna.list[[each.chr]][which(rna.list[[each.chr]]$pos %in% snp.pos.tot),'seg.b'] <- seg.pos.b
          rna.list[[each.chr]][which(rna.list[[each.chr]]$pos %in% snp.pos.tot),'seg.a'] <- seg.pos.a
        }
        if(plot.bins){
          rbin.list[[each.chr]][which(rbin.list[[each.chr]]$start.pos %in% snp.pos.tot), 'seg.a'] <- seg.pos.a
          rbin.list[[each.chr]][which(rbin.list[[each.chr]]$start.pos %in% snp.pos.tot), 'seg.b'] <- seg.pos.b
        }        
        
      }
      
    }
    chr.col.cnt <- 1 + chr.col.cnt
    dev.off()
  }
  
  dna.rna.list <- list(dna=exome.list,
                       rna=rna.list,
                       dna.bin=ebin.list,
                       rna.bin=rbin.list,
                       dna.seg=dna.seg,
                       rna.seg=rna.seg)
  return(dna.rna.list)
}


# Function: plotCbsSegments
# Purpose:  Takes a list of allelic fraction values to plot their CBS-identified segments
# Input:  chr: Current chromosome that all the SNPs belong to
#         val: allelic fraction value for every SNP given
#         pos: position of each corresponding SNP allelic fraction value
#         data.type: data type can be either 'dna' or 'rna'
#         allele: allele is either 'a' or 'b'
# Returns:  A list containing two variables - seg.pos and segmented.smoothed
#             seg.pos - vector containing the corresponding seg.mean for every SNP
#             segmented.smoothed - data.frame containing the CBS identified segments
# Plots:  Adds line segments to the allelic fraction plots for the given CBS segment
plotCbsSegments <- function(chr, val, pos, data.type, allele){
  val.pos <- c()
  seg.pos <- rep(-1, length(pos))
  
  pos.filt <- c()
  val.filt <- c()
  #pos.filt <- pos[which(val != 0 & val != 1)]
  #val.filt <- val[which(val != 0 & val != 1)]
  if(allele == 'a'){
    val.pos <- which(val <= 0.6 & val != 0 & val != 1)
    pos.filt <- pos[val.pos]
    val.filt <- val[val.pos]
  } else if (allele == 'b') {
    val.pos <- which(val >= 0.4 & val != 0 & val != 1)
    pos.filt <- pos[val.pos]
    val.filt <- val[val.pos]
  }
  
  segment.list <- list(seg.pos=seg.pos,
                       segment.smoothed=NULL)
  # Checks if there are enough values for CBS
  if(length(val.filt) >= 5){
    segment.smoothed.CNA.object <- getCbsSeg(chr, val.filt, pos.filt, paste(data.type, chr, allele, sep="-"))
    
    chr1.seg <- segment.smoothed.CNA.object$output
    for(i in seq(1:dim(chr1.seg)[1])) {
      #guideline at CNV median
      segments(chr1.seg[i,'loc.start'],
               as.numeric(chr1.seg[i,'seg.mean']),
               chr1.seg[i,'loc.end'],
               as.numeric(chr1.seg[i,'seg.mean']),
               col="orange",
               lwd=2)       
    }
    
    # Map the segments to the positions
    for(each.row in 1:dim(chr1.seg)[1]){
      row.pos <- val.pos[pos[val.pos] >= chr1.seg[each.row, 'loc.start'] & pos[val.pos] <= chr1.seg[each.row, 'loc.end']]
      seg.pos[row.pos] <- chr1.seg[each.row, 'seg.mean']
    }
    
    segment.list <- list(seg.pos=seg.pos,
                         segment.smoothed=segment.smoothed.CNA.object)
  }
  return(segment.list)
}

# Function: getCbsSeg
# Purpose:  Given the position, values, and chromosome, runs the DNAcopy CBS algorithm
# Input:  chr: Current chromosome that all the SNPs belong to
#         val: allelic fraction value for every SNP given
#         pos: position of each corresponding SNP allelic fraction value
#         name: identifier 
# Returns:  DNAcopy Object - CBS segments with $output having the data.frame
getCbsSeg <- function(chr, val, pos, name, aval=0.01, nperm=1000){
  chr.id <- unlist(regmatches(chr, regexec(pattern = "[0-9X]+", chr)))
  if(!chr.id %in% 'X'){
    chr.num <- as.numeric(chr.id)
  } else {
    chr.num <- chr.id
  }
  
  CNA.object <- CNA(val,
                    rep(chr.num, length(pos)),
                    pos, data.type="logratio", sampleid=name)
  smoothed.CNA.object <- smooth.CNA(CNA.object)
  
  segment.smoothed.CNA.object <- segment(CNA.object, verbose=1, 
                                         alpha=aval, nperm=nperm)
  
  return(segment.smoothed.CNA.object)
}


# Function: plotDensityMin
# Purpose:  Creates density plots with local min and max of densities
#     Returns the lower and upper mins for purity estimation
plotDensityMin <- function(name='undef', genomic.list, 
                           col.a, col.b, quant=0.95, min.gen.frac=0.10){
  require(scales)
  # Remove the X chromosome if present:
  if('chrX' %in% names(genomic.list)){
    genomic.list <- genomic.list[-which(names(genomic.list) %in% 'chrX')]
  }
  
  #Generate Density and remove -1 fillers with NA for removal
  genomic.df <- do.call("rbind", genomic.list)
  genomic.df[which(genomic.df$seg.a < 0), 'seg.a'] <- NA
  genomic.df[which(genomic.df$seg.b < 0), 'seg.b'] <- NA
  sega.dens <- density(genomic.df$seg.a, na.rm=TRUE)  #typically lower
  segb.dens <- density(genomic.df$seg.b, na.rm=TRUE)  #typically upper
  
  # get lower AF threshold for CBS segments
  sega.min.max <- getLocalMinMax(sega.dens, 'lower', min.gen.frac)    #get local minima for lower end of function
  lower.min <- min(sega.min.max[['min']])
  lower.af <- sega.dens$x[lower.min]
  # get upper AF threshold for CBS segments
  segb.min.max <- getLocalMinMax(segb.dens, 'upper', min.gen.frac)  #get local minima for upper end of function
  upper.min <- max(segb.min.max[['min']])
  upper.af <- segb.dens$x[upper.min]
  
  #Get local peak as purity estimate
  lower.peak <- sega.dens$x[min(sega.min.max[['max']])]
  upper.peak <- segb.dens$x[max(segb.min.max[['max']])]
  pur.est.alt <- round(getPurityEstimate(upper.peak, type='alt') * 100, 2)
  pur.est.ref <- round(getPurityEstimate(lower.peak, type='ref') * 100, 2)
  pur.est <- round((pur.est.alt + pur.est.ref)/2, 2)
  
  # Isolate raw AF for segments that meet both threshold criteria
  loh.genomic.df <- genomic.df[which(genomic.df$seg.a < lower.af & 
                                       genomic.df$seg.b > upper.af),]
  hist.a <- hist(as.numeric(as.character(loh.genomic.df[,col.a])), breaks=50, xlim=c(0,1))
  hist.b <- hist(as.numeric(as.character(loh.genomic.df[.col.b])), breaks=50, xlim=c(0,1))
  
  #Estimate the ROUGH approximation of genomic fraction attributed to each peak
  resize.factor <- 1 / sum(sega.dens$y[sega.min.max[['max']]])
  max.gen.frac <- max(sega.dens$y[sega.min.max[['max']]] * resize.factor)
  
  #rescale y density since absolute count doesn't matter
  sega.dens$y <- rescale(sega.dens$y, c(0,max.gen.frac))
  segb.dens$y <- rescale(segb.dens$y, c(0,max.gen.frac))
  
  #rescale the counts to a 0,1 scale
  hist.a$counts <- rescale(hist.a$counts, c(0,max.gen.frac))
  hist.b$counts <- rescale(hist.b$counts, c(0,max.gen.frac))
  
  #Generate CBS segment quantiles
  lower.quant <- quantile(as.numeric(as.character(loh.genomic.df$seg.a)), quant)
  upper.quant <- quantile(as.numeric(as.character(loh.genomic.df$seg.b)), (1-quant))
  
  ###### plot Densities
  pdf(file=paste(name, ".pdf", sep=""))
  split.screen(matrix(c(rep(0, 4),
                        rep(1, 4),
                        c(0.000, 0.800, 0.850, 0.900),
                        c(0.800, 0.850, 0.900, 1.000)),
                      ncol=4))
  # Adds in raw AF for SNPs included in threshold-met CBS segments
  screen(1)
  par(mar=c(5.1,4.1,0,2.1))
  plot(hist.a, col=alpha("grey", 0.4), border=alpha("grey", 0.4), 
       xlim=c(0,1), ylim=c(0, (max.gen.frac + 0.05)),
       ylab="estimated genomic fraction", xlab="segmented allelic fraction", 
       main="", sub=paste(name, ": ", pur.est, " purity", sep=""), las=1)
  plot(hist.b, col=alpha("grey", 0.4), border=alpha("grey", 0.4), add=TRUE)
  
  # Adds Density lines
  lines(sega.dens$x, sega.dens$y, xlim=c(0,1), col="black")
  lines(segb.dens$x, segb.dens$y, xlim=c(0,1), col="black")
  
  # Adds minimum genomic fraction line
  abline(h = min.gen.frac, col="black", lty=2)
  
  # Adds in minima markers
  abline(v = sega.dens$x[sega.min.max[['min']]], col="black", lty=3)
  points(sega.dens$x[sega.min.max[['min']]], sega.dens$y[sega.min.max[['min']]], col="black")
  points(sega.dens$x[lower.min], sega.dens$y[lower.min], col="black", pch=16)
  
  abline(v = segb.dens$x[segb.min.max[['min']]], col="black", lty=3)
  points(segb.dens$x[segb.min.max[['min']]], segb.dens$y[segb.min.max[['min']]],col="black")
  points(segb.dens$x[upper.min], segb.dens$y[upper.min], col="black", pch=16)
  
  # Adds in quantile thresholds that are used for purity estimation
  abline(v=c(upper.peak, lower.peak), col="red", lty=5)
  
  ###### plot Quantile Boxplots
  #Generate boxplots of CBS AF segments
  screen(2)
  par(mar=c(0,4.1,0,2.1))
  boxplot(as.numeric(as.character(loh.genomic.df$seg.a)), horizontal=TRUE, 
          xlim=c(0.75,1.25), ylim=c(0,1), xlab="test", ylab="", axes=FALSE)
  boxplot(as.numeric(as.character(loh.genomic.df$seg.b)), horizontal=TRUE, 
          xlim=c(0.75,1.25), ylim=c(0,1), axes=FALSE, add=TRUE)
  
  #Generate boxplots of raw AF segments
  screen(3)
  par(mar=c(0,4.1,0,2.1))
  boxplot(as.numeric(as.character(loh.genomic.df[,col.a])), horizontal=TRUE, 
          xlim=c(0.75,1.25), ylim=c(0,1), xlab="", ylab="", axes=FALSE)
  boxplot(as.numeric(as.character(loh.genomic.df[,col.b])), horizontal=TRUE, 
          xlim=c(0.75,1.25), ylim=c(0,1), axes=FALSE, add=TRUE)
  
  close.screen(all.screens=TRUE)
  dev.off()
  
  return(list(lowerp=lower.peak, upperp=upper.peak,
              lowerq=lower.quant, upperq=upper.quant))
}

# Function: getLocalMinMax
# Purpose:  Returns the local mininma/maxima in a kernel density function
getLocalMinMax <- function(dens.func, tail.end='upper', min.gen.frac=0.10){
  min_indexes = which(diff(  sign(diff( c(0,dens.func$y)))) == 2)
  max_indexes = which(diff(  sign(diff( c(0,dens.func$y)))) == -2)
  
  # Temporary repetitive code to identify genomic fraction according to peak height
  resize.factor <- 1 / sum(dens.func$y[max_indexes])
  max.gen.frac <- max(dens.func$y[max_indexes] * resize.factor)
  temp.y <- rescale(dens.func$y, c(0,max.gen.frac))
  
  #Removes local minima that are results of flat tails
  max_indexes <- max_indexes[temp.y[max_indexes] > min.gen.frac]
  
  if(tail.end %in% 'upper'){
    #finds the proper local minima that is lower than the maxima with the highest AF
    min_indexes <- min_indexes[dens.func$x[min_indexes] < max(dens.func$x[max_indexes])]
  } else if(tail.end %in% 'lower'){
    #finds the proper local minima that is lower than the maxima with the highest AF
    min_indexes <- min_indexes[dens.func$x[min_indexes] > min(dens.func$x[max_indexes])]
  }
  
  return(list(min=min_indexes, max=max_indexes))
}



# Function: plotDensityMin
# Purpose:  Creates density plots with local min and max of densities
#     Returns the lower and upper mins for purity estimation
plotDensityMin <- function(name='undef', genomic.list, col.a='a', col.b='b', quant=0.95, min.gen.frac=0.10){
  require(scales)
  #Generate Density and remove -1 fillers with NA for removal
  genomic.df <- do.call("rbind", genomic.list)
  genomic.df[which(genomic.df$seg.a < 0), 'seg.a'] <- NA
  genomic.df[which(genomic.df$seg.b < 0), 'seg.b'] <- NA
  sega.dens <- density(genomic.df$seg.a, na.rm=TRUE)  #typically lower
  segb.dens <- density(genomic.df$seg.b, na.rm=TRUE)  #typically upper
  
  # get lower AF threshold for CBS segments
  sega.min.max <- getLocalMinMax(sega.dens, 'lower', min.gen.frac)    #get local minima for lower end of function
  lower.min <- min(sega.min.max[['min']])
  lower.af <- sega.dens$x[lower.min]
  # get upper AF threshold for CBS segments
  segb.min.max <- getLocalMinMax(segb.dens, 'upper', min.gen.frac)  #get local minima for upper end of function
  upper.min <- max(segb.min.max[['min']])
  upper.af <- segb.dens$x[upper.min]
  
  #Get local peak as purity estimate
  lower.peak <- sega.dens$x[min(sega.min.max[['max']])]
  upper.peak <- segb.dens$x[max(segb.min.max[['max']])]
  pur.est.alt <- round(getPurityEstimate(upper.peak, type='alt') * 100, 2)
  pur.est.ref <- round(getPurityEstimate(lower.peak, type='ref') * 100, 2)
  pur.est <- round((pur.est.alt + pur.est.ref)/2, 2)
  
  # Isolate raw AF for segments that meet both threshold criteria
  loh.genomic.df <- genomic.df[which(genomic.df$seg.a < lower.af & 
                                       genomic.df$seg.b > upper.af),]
  hist.a <- hist(as.numeric(as.character(loh.genomic.df[,col.a])), breaks=50, xlim=c(0,1))
  hist.b <- hist(as.numeric(as.character(loh.genomic.df[,col.b])), breaks=50, xlim=c(0,1))
  
  #Estimate the ROUGH approximation of genomic fraction attributed to each peak
  resize.factor <- 1 / sum(sega.dens$y[sega.min.max[['max']]])
  max.gen.frac <- max(sega.dens$y[sega.min.max[['max']]] * resize.factor)
  
  #rescale y density since absolute count doesn't matter
  sega.dens$y <- rescale(sega.dens$y, c(0,max.gen.frac))
  segb.dens$y <- rescale(segb.dens$y, c(0,max.gen.frac))
  
  #rescale the counts to a 0,1 scale
  hist.a$counts <- rescale(hist.a$counts, c(0,max.gen.frac))
  hist.b$counts <- rescale(hist.b$counts, c(0,max.gen.frac))
  
  #Generate CBS segment quantiles
  lower.quant <- quantile(as.numeric(as.character(loh.genomic.df$seg.a)), quant)
  upper.quant <- quantile(as.numeric(as.character(loh.genomic.df$seg.b)), (1-quant))
  
  ###### plot Densities
  pdf(file=paste(name, ".pdf", sep=""))
  split.screen(matrix(c(rep(0, 4),
                        rep(1, 4),
                        c(0.000, 0.800, 0.850, 0.900),
                        c(0.800, 0.850, 0.900, 1.000)),
                      ncol=4))
  # Adds in raw AF for SNPs included in threshold-met CBS segments
  screen(1)
  par(mar=c(5.1,4.1,0,2.1))
  plot(hist.a, col=alpha("grey", 0.4), border=alpha("grey", 0.4), 
       xlim=c(0,1), ylim=c(0, (max.gen.frac + 0.05)),
       ylab="estimated genomic fraction", xlab="segmented allelic fraction", 
       main="", sub=paste(name, ": ", pur.est, " purity", sep=""), las=1)
  plot(hist.b, col=alpha("grey", 0.4), border=alpha("grey", 0.4), add=TRUE)
  
  # Adds Density lines
  lines(sega.dens$x, sega.dens$y, xlim=c(0,1), col="black")
  lines(segb.dens$x, segb.dens$y, xlim=c(0,1), col="black")
  
  # Adds minimum genomic fraction line
  abline(h = min.gen.frac, col="black", lty=2)
  
  # Adds in minima markers
  abline(v = sega.dens$x[sega.min.max[['min']]], col="black", lty=3)
  points(sega.dens$x[sega.min.max[['min']]], sega.dens$y[sega.min.max[['min']]], col="black")
  points(sega.dens$x[lower.min], sega.dens$y[lower.min], col="black", pch=16)
  
  abline(v = segb.dens$x[segb.min.max[['min']]], col="black", lty=3)
  points(segb.dens$x[segb.min.max[['min']]], segb.dens$y[segb.min.max[['min']]],col="black")
  points(segb.dens$x[upper.min], segb.dens$y[upper.min], col="black", pch=16)
  
  # Adds in quantile thresholds that are used for purity estimation
  abline(v=c(upper.peak, lower.peak), col="red", lty=5)
  
  ###### plot Quantile Boxplots
  #Generate boxplots of CBS AF segments
  screen(2)
  par(mar=c(0,4.1,0,2.1))
  boxplot(as.numeric(as.character(loh.genomic.df$seg.a)), horizontal=TRUE, 
          xlim=c(0.75,1.25), ylim=c(0,1), xlab="test", ylab="", axes=FALSE)
  boxplot(as.numeric(as.character(loh.genomic.df$seg.b)), horizontal=TRUE, 
          xlim=c(0.75,1.25), ylim=c(0,1), axes=FALSE, add=TRUE)
  
  #Generate boxplots of raw AF segments
  screen(3)
  par(mar=c(0,4.1,0,2.1))
  boxplot(as.numeric(as.character(loh.genomic.df[,col.a])), horizontal=TRUE, 
          xlim=c(0.75,1.25), ylim=c(0,1), xlab="", ylab="", axes=FALSE)
  boxplot(as.numeric(as.character(loh.genomic.df[,col.b])), horizontal=TRUE, 
          xlim=c(0.75,1.25), ylim=c(0,1), axes=FALSE, add=TRUE)
  
  close.screen(all.screens=TRUE)
  dev.off()
  
  return(list(lowerp=lower.peak, upperp=upper.peak,
              lowerq=lower.quant, upperq=upper.quant))
}


# Function: getPurityEstimate
# Purpose:  Calculates a purity estimate given an alternate allele AF
getPurityEstimate <- function(af.peak=NULL, type='alt', pl.t=2, f.t=1, 
                              pl.n=2, f.n=0.5){
  if(type %in% 'alt'){
    pur.t <- ((pl.n*(f.n - af.peak)) / 
                (pl.t*f.t*(af.peak-1) + pl.n*(f.n-af.peak)))  
  } else if (type %in% 'ref'){
    pur.t <- ((pl.n*(f.n + af.peak - 1))/
                ((pl.n*(f.n + af.peak - 1)) - (pl.t*f.t*af.peak)))
  }
  
  return(pur.t)
}

# Function: plotTNnormalization
# Purpose: Plots the Tumor - Normal normalization difference in bin snp density
plotTNnormalization <- function(seg.x.df, seg.y.df){
  pdf("segDensNormalization.pdf")
  split.screen(c(1,2))
  screen(1)
  plot(density(seg.x.df$seg.ab, na.rm=TRUE), xlab="allelic fraction", main="raw AF density")
  lines(x=density(seg.y.df$seg.ab, na.rm=TRUE)$x,
        y=density(seg.y.df$seg.ab, na.rm=TRUE)$y, col="red")  
  legend("topright", legend=c("tumor", "normal"), col=c("black", "red"), lty=c(1,1), cex=0.8)

  screen(2)
  plot(density(seg.x.df$seg.ab.norm, na.rm=TRUE), xlab="allelic fraction", main="median normalized AF density")
  lines(x=density(seg.y.df$seg.ab.norm, na.rm=TRUE)$x,
        y=density(seg.y.df$seg.ab.norm, na.rm=TRUE)$y, col="red")
  legend("topright", legend=c("tumor", "normal"), col=c("black", "red"), lty=c(1,1), cex=0.8)
  close.screen(all.screens=TRUE)
  dev.off()  
}

# Function: plotTumorNormDiff
# Purpose: Plots the mean segment value and associated standard deviation
#   between tumor and normal samples.  As well as the associated bonferroni
#   corrected p value
plotTumorNormDiff <- function(tn.df, suffix=""){
  tn.df$chr <- as.character(tn.df$chr)
  tn.df$chr <- factor(tn.df$chr, levels=paste("chr", c(1:22), sep=""))
  tn.split <- split(tn.df, tn.df$chr)
 
  pdf(paste("T-N.BinDifference", suffix, ".pdf", sep=""), height=30)
  chr.screens <- split.screen(c(length(tn.split), 1))
  lapply(chr.screens, function(chr.cnt){
    screen(chr.cnt)
    x <- tn.split[[chr.cnt]]
    
    arm.split <- split(x, x$arm)
    arm.screens <- split.screen(c(1,2))
    if(length(arm.split) == 1){
      arm.split <- append(list(arm.split[[1]][1,]), arm.split)
      arm.split[[1]][1, c('mean.seg.T', 'mean.sd.T', 'mean.seg.N', 'mean.sd.N')] <- NA
    }
    for(each.arm in c(1,2)){
      # Format the arm-level dataframe
      arm.df <- arm.split[[each.arm]]
      arm.df <- arm.df[which(!is.na(arm.df$pval)),]
      
      screen(arm.screens[each.arm])
      if(each.arm==1){
        par(mar=c(0.1, 4.1, 0.1, 0.1))
      } else if (each.arm == 2) {
        par(mar=c(0.1, 0, 0.1, 4.1))
      }
      plot(0, type='n', xlim=c(arm.split[[each.arm]][1,'start.pos'],
                               arm.split[[each.arm]][dim(arm.split[[each.arm]])[1],'end.pos']), 
           ylim=c(0,0.5), yaxt=if(each.arm==1) 't' else 'n', xaxt='n', 
           ylab=if(each.arm==1) names(tn.split)[chr.cnt] else "", xlab="", las=2, cex=0.8)
      
      rect(xleft=arm.df$start.pos, ybottom=(arm.df$mean.seg.T - arm.df$mean.sd.T),
           xright=arm.df$end.pos, ytop=(arm.df$mean.seg.T + arm.df$mean.sd.T), 
           col=alpha("cyan", 0.1), lty=0)
      rect(xleft=arm.df$start.pos, ybottom=(arm.df$mean.seg.N - arm.df$mean.sd.N),
           xright=arm.df$end.pos, ytop=(arm.df$mean.seg.N + arm.df$mean.sd.N), 
           col=alpha("cyan", 0.1), lty=0)
      
      
#       pvals <- arm.df$pval.adj
#       pvals[which(pvals > 0.2)] <- 0.2
#       pval.index <- as.integer(rescale(pvals, to=c(0,1000)))
#       pvalColorT <- colorRampPalette(c("Orange", "Orange", "grey"))(1000)
#       pvalColorN <- colorRampPalette(c("Black", "Black", "grey"))(1000)
      pval.index <- which(arm.df$pval.adj < 0.05)
      pvalColorT <- rep("grey", dim(arm.df)[1])
      pvalColorT[pval.index] <- rep("Orange", length(pval.index))
      pvalColorN <- gsub("Orange", "Black", pvalColorT)
      rect(xleft=arm.df$start.pos, ybottom=(arm.df$mean.seg.T - 0.01), 
           xright=arm.df$end.pos, ytop=(arm.df$mean.seg.T  + 0.01), 
           col=pvalColorT, lty=0)
      rect(xleft=arm.df$start.pos, ybottom=(arm.df$mean.seg.N - 0.01), 
           xright=arm.df$end.pos, ytop=(arm.df$mean.seg.N  + 0.01), 
           col=pvalColorN, lty=0)
    }
  })
  close.screen(all.screens=TRUE)
  dev.off()
}

# Function: getSigmaEstimate
# Purpose: With a list of density maximas and minimas, it finds the minima corresponding to each
#   peak maxima.  Then takes the average of the distance between the minimas and maxima
#   and takes that as representative of 3 stdDev (99.5% quantile).  Divides by 3 to get sigma estimate
getSigmaEstimate <- function(lmax, lmin, max.val=0.5, min.val=-1000){
  require(intervals)
  if(min(lmin) > min(lmax)){
    lmin <- c(min.val, lmin, max.val)
  } else {
    lmin <- c(lmin, max.val)
  }
  
  # Create intervals of start and end points for each peak based on minimas, and the maxima point
  min.mat <- matrix(c(lmin[-length(lmin)], lmin[-1]), byrow=FALSE, ncol=2)
  min.int <- Intervals(min.mat, closed=TRUE)
  max.int <- Intervals(matrix(c(lmax, lmax), byrow=FALSE, ncol=2), closed=TRUE)
  
  # Create a peak matrix with the local-minimas corresponding to each maxima
  peak.mat <- cbind(lmax, min.mat[unlist(interval_overlap(max.int, min.int)),])
  
  # Get the sigma estimates
  peak.sigma <- apply(peak.mat, 1, function(x){
    peak.disp <- (abs(x[2] - x[1]) + abs(x[3] - x[1])) / 2
    peak.sigma <- peak.disp / 3   # Assumes the dispersion limits are representative of 3 stdDev
    return(peak.sigma)
  })
  return(peak.sigma)
}

# Function: decomposeAfDist
# Purpose:  Uses the MIX algorithm to apply EM and Newton-methods to decompose
#   a mixed multi-gaussian allelic fraction profile
decomposeAfDist <- function(genomic.list,  col.a, col.b){
  require(scales)
  require(mixdist)
  # Remove the X chromosome if present:
  if('chrX' %in% names(genomic.list)){
    genomic.list <- genomic.list[-which(names(genomic.list) %in% 'chrX')]
  }
  
  #Generate Density and remove -1 fillers with NA for removal
  genomic.df <- do.call("rbind", genomic.list)
  genomic.df[which(genomic.df$seg.a < 0), 'seg.a'] <- NA
  genomic.df[which(genomic.df$seg.b < 0), 'seg.b'] <- NA
  seg.tot <- c(abs(genomic.df$seg.a), abs(1 - genomic.df$seg.b))
  seg.tot.dens <- density(seg.tot, na.rm=TRUE)  #typically lower
  
  # get lower AF threshold for CBS segments
  seg.min.max <- getLocalMinMax(seg.tot.dens, 'lower', 0.05)    #get local minima for lower end of function
  
  break.num <- 30
  seg.his <-  hist(seg.tot, breaks=break.num)
  segAf.df <- data.frame(mid=seg.his$mids, cou=seg.his$counts)  
  
  
  estMean <- seg.tot.dens$x[seg.min.max$max]
  estSigma <- getSigmaEstimate(estMean, seg.tot.dens$x[seg.min.max$min])
  estDist <- 'norm'
  segAf.df <- segAf.df[which(segAf.df$mid < (estMean[length(estMean)] + (4 * estSigma[length(estSigma)]))),]
  
  fitpro <- mix(as.mixdata(segAf.df), mixparam(mu=estMean, sigma=estSigma), dist=estDist)
  return(fitpro)
}
  
