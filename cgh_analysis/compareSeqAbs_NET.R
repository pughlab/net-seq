####################################
###### aggregateLoh
######   Purpose: To cgroup together LOH profiles and 
######      total copy-number information between samples
######      of the same disease type (i.e. pnets).
######   Usage:  Takes the output of plotRawLoH.R as the
######      input.
######         net.list <- list(pnets=pnet.seg.list, 
######                          sinets=sinet.seg.list,
######                          ccl=pnet.ccl.seg.list)
######         num.of.samples <- sum(unlist(lapply(net.list, function(x) length(names(x)))))
######         disease.list <- net.list
####################################
library(intervals)
library(scales)

###############################
#         Variables
###############################
remove.sex <- 1   # Removes the chr X and Y
ucsc.chrom.info.file <- '/Users/rquevedo/Desktop/bhk_lab/reference/ucsc.hg19.chromInfo.txt'     #hg19

#Sample Info
pnet.seg.file <- '/Users/rquevedo/Desktop/PughLab/NET-seq/Sequenza/bkup/all/segments/PNET_segments.txt'
sample.id <- 'pancreatic NET - Sequenza'
abs.pnet.seg.file <- '/Users/rquevedo/Desktop/PughLab/NET-seq/snp6/absolute/output/absolute/abs_extract/reviewed/SEG_MAF/NET.all.segtab.txt'
sample.id2 <- 'pancreatic NET - Absolute'

disease.names <- c("net001", "net003", "net008", "net009")
each.disease <- disease.names[1]
dis.cn.pal <- c('white', 'dodgerblue', 'thistle','red') # CN = 0, CN = 1, CN = 2, CN = 3+ 

###############################
#         Functions
###############################
# Matches a given samples LOH segment to a master interval breakpoint list to retrieve the respective total copy-number
getCNt <- function(sample.int.df, st.bp, end.bp, cn.t, each.dis){
  single.seg.int <- Intervals(matrix(c(st.bp, end.bp), ncol=2), closed=TRUE)
  if(cn.t > 3) cn.t <- 3
  sample.int.df[which(interval_overlap(t.int, single.seg.int) != 0),each.dis] <- cn.t
  sample.int.df <- sample.int.df[which(sample.int.df[,each.dis] != 0),]
  return(sample.int.df)
}

getChrBp <- function(ucsc.chrom, x){
  x <- as.integer(as.character(x))
  chr.lengths <- ucsc.chrom[which(ucsc.chrom$chr.num < x), 'size']
  chr.size <- sum(as.numeric(as.character(chr.lengths)))
  return(chr.size)
}



###############################
#           Main
###############################
########## Pre Process the Data ##########
#Read in and obtain NET Segments
pnet.seg <- read.csv(pnet.seg.file, header = TRUE, sep = "\t", quote = "\"", check.names=FALSE, stringsAsFactors=FALSE)
pnet.seg.list <- split(pnet.seg, pnet.seg$Sample)
pnet.seg.list[['Sample']] <- NULL

abs.pnet.seg <- read.csv(abs.pnet.seg.file, header = TRUE, sep = "\t", quote = "\"", check.names=FALSE, stringsAsFactors=FALSE)
abs.pnet.seg$Chromosome <- gsub("^", "chr", abs.pnet.seg$Chromosome)
colnames(abs.pnet.seg)[which(colnames(abs.pnet.seg) %in% c("sample", "modal.a1", "modal.a2"))] <- c("Sample", "modal_A1", "modal_A2")
abs.pnet.seg$CNt <- (as.integer(abs.pnet.seg$modal_A1) + as.integer(abs.pnet.seg$modal_A2))
abs.pnet.seg.list <- split(abs.pnet.seg, abs.pnet.seg$Sample)

# Manually adjust NET-008 ploidy calls
abs.net <- abs.pnet.seg.list[['NET008FF2']]
abs.net$modal_A1 <- (as.integer(abs.net$modal_A1) - 1)
abs.net$modal_A2 <- (as.integer(abs.net$modal_A2) - 1)
abs.net$CNt <- (abs.net$CNt - 1)
abs.net[which(abs.net$modal_A1 %in% -1), 'modal_A1'] <- 0
abs.net[which(abs.net$modal_A1 %in% 0), 'LOH'] <- 1
abs.pnet.seg.list[['NET008FF2']] <- abs.net

# Manually adjust NET-009 ploidy calls
abs.net <- abs.pnet.seg.list[['NET009FF1']]
abs.net$modal_A1 <- (as.integer(abs.net$modal_A1) - 1)
abs.net$modal_A2 <- (as.integer(abs.net$modal_A2) - 1)
abs.net$CNt <- (abs.net$CNt - 1)
abs.net[which(abs.net$modal_A1 %in% -1), 'modal_A1'] <- 0
abs.net[which(abs.net$modal_A1 %in% 0), 'LOH'] <- 1
abs.pnet.seg.list[['NET009FF1']] <- abs.net

# Remove the sample from list that is an artifact from poor formatting. sorry life.
abs.pnet.seg.list[['sample']] <- NULL

# Create a master list for seq-abs paired segments
pnet.list <- list(net001=list(net001.seq=pnet.seg.list[['NET-001-T_2']],
                              net001.abs=abs.pnet.seg.list[['NET001FF1']]), 
                  net003=list(net003.seq=pnet.seg.list[['NET-003-T_2']],
                              net003.abs=abs.pnet.seg.list[['NET003FF3']]), 
                  net008=list(net008.seq=pnet.seg.list[['NET-008-T_2']],
                              net008.abs=abs.pnet.seg.list[['NET008FF2']]), 
                  net009=list(net009.seq=pnet.seg.list[['NET-009-T_2']],
                              net009.abs=abs.pnet.seg.list[['NET009FF1']]))
num.of.samples <- sum(unlist(lapply(pnet.list, function(x) length(names(x)))))
disease.list <- pnet.list

# read and filter ucsc.chrom.info file
ucsc.chrom <- read.csv(ucsc.chrom.info.file, header=TRUE, sep="\t", check.names=FALSE)
colnames(ucsc.chrom) <- gsub("#", "", colnames(ucsc.chrom))
ucsc.chrom <-  ucsc.chrom[regexpr("_", ucsc.chrom$chrom) == -1,]
ucsc.chrom <- ucsc.chrom[regexpr("^chrm$", ucsc.chrom$chrom, ignore.case = TRUE, perl=TRUE) == -1,]
if(remove.sex == 1){
  ucsc.chrom <- ucsc.chrom[regexpr("^chr[xy]$", ucsc.chrom$chrom, ignore.case = TRUE, perl=TRUE) == -1,]
}
#Calculate the chromosomal pseudo start bp for plots
chr.num <- as.numeric(as.character(gsub("^chr", "", ucsc.chrom$chrom)))
ucsc.chrom <- cbind(ucsc.chrom, chr.num=chr.num)
apply(ucsc.chrom, 1, function(x) getChrBp(ucsc.chrom, x['chr.num']))
ucsc.chrom <- cbind(ucsc.chrom, chr.st.size = apply(ucsc.chrom, 1, function(x) getChrBp(ucsc.chrom, x['chr.num'])))
total.genome.size <- sum(as.numeric(ucsc.chrom$size))




########## VISUALIZATION ##########
# Creates a legend for addition to the graphs
print.legend <- 0
if(print.legend == 1){
  pdf("legends.pdf")
  # Prints a legend
  plot (0, type='n', xlim=c(0,10), ylim=c(0,10), axes=FALSE, ylab="", xlab="")
  legend(1,9, 
         c("Copy-Loss LOH","Copy-Neutral LOH", "Copy-Gain LOH"),
         inset=0.02,
         border="black",
         fill=dis.cn.pal[c(2:4)],
         cex=0.5) 
  dev.off()
  
}

setwd('/Users/rquevedo/Desktop/')
pdf("chr.net.all.pdf", width=13)
# pnets - seq:abs layout
layout(matrix(c(rep(0,30),rep(0,30),
                0,45,45,45,seq(from=1, to=44, by=2),rep(0,4),rep(0,30),
                0,45,45,45,seq(from=1, to=44, by=2),rep(0,4),rep(0,30),
                0,46,46,46,seq(from=2, to=44, by=2),rep(0,4),rep(0,30),
                rep(0,30),rep(0,30),
                0,91,91,91,seq(from=47, to=90, by=2),rep(0,4),rep(0,30),
                0,91,91,91,seq(from=47, to=90, by=2),rep(0,4),rep(0,30),
                0,92,92,92,seq(from=48, to=90, by=2),rep(0,4),rep(0,30),
                rep(0,30),rep(0,30),
                0,137,137,137,seq(from=93, to=136, by=2),rep(0,4),rep(0,30),
                0,137,137,137,seq(from=93, to=136, by=2),rep(0,4),rep(0,30),
                0,138,138,138,seq(from=94, to=136, by=2),rep(0,4),rep(0,30),
                rep(0,30),rep(0,30),
                0,183,183,183,seq(from=139, to=182, by=2),rep(0,4),rep(0,30),
                0,183,183,183,seq(from=139, to=182, by=2),rep(0,4),rep(0,30),
                0,184,184,184,seq(from=140, to=182, by=2),rep(0,4),rep(0,30)), 
              ncol=16, byrow=FALSE))

disease.names <- c("net001", "net003", "net008", "net009")
for(each.disease in disease.names){
  all.chr.int.df <- matrix(ncol=(4 + length(names(disease.list[[each.disease]]))), nrow=0)
  colnames(all.chr.int.df) <- c("chrom", "Start.bp", "End.bp", names(disease.list[[each.disease]]), "t_loh")
  
  for(each.chr in gsub("^", "chr", c(1:22))){
    print(paste("Dealing with ", each.disease, ":", each.chr, sep=""))
    
    ### Goes through each sample and creates an Interval list for each chromosome
    chr.matrix.list <- list() # Master Intervals list
    for(each.sample in names(disease.list[[each.disease]])){
      sample.df <- disease.list[[each.disease]][[each.sample]]
      
      ####################### TEMPORARY FIX #######################
      # Make sure a total CN column is present - CNt
      sample.df$CNt <- as.integer(as.character(sample.df$modal_A1)) + as.integer(as.character(sample.df$modal_A2))
      if(length(grep("^[0-9]", sample.df$Chromosome)) > 0){
        sample.df$Chromosome <- gsub("^", "chr", sample.df$Chromosome)
      }
      ####################### \TEMPORARY FIX ######################
      
      # formatting class & subsetting
      sample.df$Start.bp <- as.integer(as.character(sample.df$Start.bp))
      sample.df$End.bp <- as.integer(as.character(sample.df$End.bp))
      
      sample.df <- sample.df[which(each.chr == sample.df$Chromosome),]
      sample.df <- sample.df[which(sample.df$LOH == 1),]
      
      #Appending intervals to the master Intervals list
      if(dim(sample.df)[1] > 0){
        sample.chr.int <- Intervals(as.matrix(sample.df[,c("Start.bp", "End.bp")]), closed=TRUE)
        
        chr.matrix.list[[each.sample]] <- sample.chr.int  
      } else {
        chr.matrix.list[[each.sample]] <- Intervals() #append a blank interval if no LOH found on chromosome
      }
    }
    
    ucsc.chr.row <- ucsc.chrom[which(ucsc.chrom$chrom %in% each.chr),]
    
    ### Create's a set of Intervals that accounts for all Samples LOH Intervals to allow for tallying of occurences
    # Create a list of ALL intersections between all sample Intervals
    total.matrix.int <- matrix(ncol=2)
    for(each.l in chr.matrix.list){
      for(pre.l in chr.matrix.list){
        total.matrix.int <- rbind(total.matrix.int, as.matrix(interval_intersection(each.l, pre.l)))
      }
    }
    # Create a non-repetitive interval list between ALL intersection points
    total.breakpoints <- unique(sort(total.matrix.int))
    total.int <- matrix(c(total.breakpoints[-length(total.breakpoints)],
                          total.breakpoints[-1]), ncol=2, byrow=F)
    t.int <- Intervals(total.int, closed=c(FALSE,FALSE))  # t.int = master intervals list for all intersections
    
    # Attribute each of these intervals to the original intervals for each disease
    total.int.df <- as.data.frame(total.int)
    colnames(total.int.df) <- c("Start.bp", "End.bp")
    if(dim(total.int.df)[1] > 0){
      for(each.dis in names(chr.matrix.list)){
        sample.df <- disease.list[[each.disease]][[each.dis]]
        ####################### TEMPORARY FIX #######################
        # Make sure a total CN column is present - CNt
        sample.df$CNt <- as.integer(as.character(sample.df$modal_A1)) + as.integer(as.character(sample.df$modal_A2))
        if(length(grep("^[0-9]", sample.df$Chromosome)) > 0){
          sample.df$Chromosome <- gsub("^", "chr", sample.df$Chromosome)
        }
        ####################### \TEMPORARY FIX ######################
        sample.df <- sample.df[which(each.chr == sample.df$Chromosome),]
        sample.df <- sample.df[which(sample.df$LOH == 1),]
        
        if(dim(sample.df)[1] > 0){
          total.int.df[,each.dis] <- 0
          # Create a list of all CN.total values for each segment
          cn.sample.df <- apply(sample.df, 1, function(x) getCNt(total.int.df,
                                                                 as.integer(as.character(x['Start.bp'])), 
                                                                 as.integer(as.character(x['End.bp'])),
                                                                 as.integer(as.character(x['CNt'])),
                                                                 each.dis))
          cn.sample.df <- do.call("rbind", cn.sample.df)
          #Replace the 0's with the proper CN.t (1,2,3 [where 3 is 3+])
          total.int.df[which(total.int.df$Start.bp %in% cn.sample.df$Start.bp),each.dis] <- cn.sample.df[,each.dis]
        } else {
          total.int.df[,each.dis] <- 0
        }
      }
    } else {
      total.int.df <- as.data.frame(matrix(ncol=(dim(total.int.df)[2] + length(names(chr.matrix.list))),
                                           nrow=0))
      colnames(total.int.df) <- c("Start.bp", "End.bp", names(chr.matrix.list))
    }
    # Creates a tally of occurence in the samples
    total.int.df[,'t_loh'] <- apply(total.int.df, 1, function(x) length(which(x[c(names(chr.matrix.list))] > 0)))
    total.int.df <- cbind(chrom=rep(each.chr, dim(total.int.df)[1]), total.int.df)
    
    num.of.dis <- length(names(chr.matrix.list))
    #  PLOT 1 - Plots the LOH fragment for EACH chromosome and EACH sample
    par(mar=c(0,4.1,0,0))
    plot(0, type='n', 
         xlim=c(0.5,(num.of.dis + 0.5)), ylim=c(1,ucsc.chr.row$size), 
         xaxt='n', yaxt='n',
         ylab='', xlab='')
    #Create a function to generate a continuous color palette
    dis.cnt <- 1
    for(each.dis in names(chr.matrix.list)){  
      each.dis.df <- total.int.df[which(total.int.df[,each.dis] >= 1),]
      if(dim(each.dis.df)[1] > 0){
        dis.cn.col <- each.dis.df[,each.dis]
        for(each.cn.t in c(0:3)){
          dis.cn.col <- gsub(paste("^", each.cn.t, "$", sep=""), dis.cn.pal[each.cn.t + 1], dis.cn.col) 
        }
        
        rect(xleft = (dis.cnt - 0.5), ybottom = each.dis.df$Start.bp,
             xright = (dis.cnt + 0.5), ytop = each.dis.df$End.bp,
             col=dis.cn.col, border=NA)  
      }
      
      dis.cnt <- dis.cnt + 1
    }
    abline(v = seq(from=0.5, to=(dis.cnt - 0.5), by = 1))
    axis(side = 2, at=(ucsc.chr.row$size / 2), labels = gsub("^chr", "", each.chr), las=2)
    if(each.chr == 'chr22'){
      axis(side = 1, at = c(1:(dis.cnt - 1)), labels = names(chr.matrix.list), las=2)
    }
    
    
    #  PLOT 2 - Plots the LOH fragment representing the aggregation of all samples
    par(mar=c(0,0,0,2.1))
    plot(0, type='n', 
         xlim=c(1,1), ylim=c(1,ucsc.chr.row$size),
         xaxt='n', yaxt='n',
         ylab='', xlab='')
    #Create a function to generate a continuous color palette
    disPal <- colorRampPalette(c('white','gray22'))
    disPal <- disPal(num.of.dis + 1)
    dis.col <- total.int.df$t_loh
    for(each.loh.cnt in c(0:num.of.dis)){
      dis.col <- gsub(each.loh.cnt, disPal[each.loh.cnt + 1], dis.col) 
    }
    
    
    if(dim(total.int.df)[1] > 0){
      rect(xleft = 0.5, ybottom = total.int.df$Start.bp,
           xright = 1.5, ytop = total.int.df$End.bp,
           col=dis.col, border=NA)
    }
    
    if(each.chr == 'chr22'){
      axis(side = 1, at = 1, labels = 'Total LOH', las=2)
    }
    
    all.chr.int.df <- rbind(all.chr.int.df, total.int.df)
  }
  
  # Calculates the percent of the total genome that is LOH
  all.chr.int.df[,'length'] <- all.chr.int.df[,'End.bp'] - all.chr.int.df[,'Start.bp']
  all.frac.loh <- c()
  for(each.sample.name in names(disease.list[[each.disease]])){
    print(each.sample.name)
    sample.df <- all.chr.int.df[which(all.chr.int.df[,each.sample.name] > 0),]
    sample.frac.loh <-  ( sum(sample.df$length) / total.genome.size ) * 100
    all.frac.loh <- c(all.frac.loh, sample.frac.loh) 
  }
  names(all.frac.loh) <- names(disease.list[[each.disease]])
  
  par(mar=c(0.25,4.1,0,0))
  barplot(all.frac.loh, space=0, xlim=c(0,length(names(disease.list[[each.disease]]))), 
          ylim=c(0,100), las=2, ylab="% LOH", xaxt='n')
  
  # Barplot for the total_loh for all samples to none
  t.loh.samples <- do.call("rbind", lapply(split(all.chr.int.df, all.chr.int.df$t_loh), function(x) (sum(x$length) / total.genome.size) * 100))
  t.loh.samples[1,1] <- (100 - sum(t.loh.samples[c(2:dim(t.loh.samples)[1]),]))
  t.loh.samples <- apply(t.loh.samples, 2, rev)
  
  disPal <- colorRampPalette(c('white','gray22'))
  disPal <- disPal(num.of.dis + 1)
  dis.col <- rownames(t.loh.samples)
  for(each.loh.cnt in c(0:num.of.dis)){
    dis.col <- gsub(each.loh.cnt, disPal[each.loh.cnt + 1], dis.col) 
  }
  
  par(mar=c(0.25,0,0,2.0))
  barplot(t.loh.samples, ylim=c(0,100), col=dis.col, yaxt='n', xaxt='n')
  
  print(paste("Disease: ", each.disease, sep=""))
  print(all.frac.loh)
  print(t.loh.samples)


  # Output the LOH segment information into .seg files for IGV plotting
  out.seg.df <- all.chr.int.df
  out.seg.df$chrom <- gsub('^chr', '', out.seg.df$chrom)
  out.seg.df$ID <- rep(paste(each.disease, "-absolute", sep=""), dim(out.seg.df)[1])
  colnames(out.seg.df) <- c("chrom", "loc.start", "loc.end", "seg.mean.sequenza", "seg.mean.absolute", "seg.mean.intersect", "num.mark", "ID")
  for(platform in c("abs", "seq", "int")){
    if(platform %in% 'abs'){
      plt <- 'seg.mean.absolute'
      id <- paste(each.disease, "-affy6", sep="")
    } else if (platform %in% 'seq'){
      plt <- 'seg.mean.sequenza'
      id <- paste(each.disease, "-wes", sep="")
    } else if (platform %in% "int"){
      plt <- 'seg.mean.intersect'
      id <- paste(each.disease, "-intersect", sep="")
    }
    plat.seg.df <- out.seg.df[,c('ID', 'chrom', 'loc.start', 'loc.end', 'num.mark', plt)]
    plat.seg.df$ID <- id
    if(platform %in% 'int'){
      plat.seg.df <- plat.seg.df[which(plat.seg.df[,plt] == 2),]  
    } else {
      plat.seg.df <- plat.seg.df[which(plat.seg.df[,plt] != 0),]  
    }
    plat.seg.df[,plt] <- 1
    write.table(plat.seg.df, file=paste(each.disease, platform, "seg", sep="."), col.names=TRUE, row.names=FALSE, quote=FALSE, sep="\t")
  }
  
}
dev.off()
#save(disease.list, total.genome.size, ucsc.chrom, file="data.aggregateLoh.Rdata")
