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

###############################
#         Variables
###############################
#save(ucsc.chrom, total.genome.size, file="/Users/rquevedo/git/loh_cn_visualize/data/chromInfo.Rdata")
#load("/Users/rquevedo/git/shallowWgsCnv/data/otbData.aggregateLoh.medianLOH.Rdata")  # otb-PNET, otb-GINET
#load("/Users/rquevedo/git/shallowWgsCnv/data/otbData.aggregateLoh.muLOH.Rdata")  # otb-PNET, otb-GINET

load("/Users/rquevedo/git/loh_cn_visualize/data/chromInfo.Rdata")
#load("/Users/rquevedo/Desktop/net_OTB/plots/otbData.aggregateLoh.medianLOH.Rdata")  # otb-PNET, otb-GINET
#load("/Users/rquevedo/git/loh_cn_visualize/data/data.aggregateLoh.Rdata")  #PNETs, GINETs, CCL   (T_2 only)
#load("~/git/loh_cn_visualize/data/data.aggregateLoh.T2-1.Rdata") #PNETs, GINETs, CCL  (T_2 and T_1)

disease.names <- c("pnets", "sinets", "ccl")
disease.names <- c("otb-pnet", "otb-ginet")
total.disease.chr.list <- list()
each.disease <- disease.names[2]
dis.cn.pal <- c('white', 'dodgerblue', 'thistle','red') # CN = 0, CN = 1, CN = 2, CN = 3+ 
mk.plots <- TRUE

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



###############################
#           Main
###############################
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
pdf(paste("chr.", each.disease, ".pdf", sep=""))


##  Setup the layout
if(any(each.disease %in% c('pnets', 'sinets', 'ccl'))){
  print(paste(each.disease, " layout", sep=""))
  padding.col <- 2
  num.samples <- length(disease.list[[each.disease]])
} else if(any(each.disease %in% c('otb-pnet', 'otb-ginet'))){
  print(paste(each.disease, " layout", sep=""))
  padding.col <- 2
  num.samples <- length(disease.list[[each.disease]])
  #order the diseases in ascending order
  if(mk.plots){
    if(each.disease %in% 'otb-pnet'){
      #order(all.frac.loh)
      #med 0.20 c(12,7,11,10,1,6,9,2,4,3,5,8)
      #med.seg 0.25 c(10,4,3,14,11,12,9,2,6,8,7,5,1,13)
      #med.seg.stroma c(9,3,10,11,13,8,2,5,7,1,4,6,12)
      # mu c(3,4,10,11,14,12,9,2,13,6,7,8,5,1)
      # mu.seg.0.25 c(10,3,4,9,14,13,12,11,5,6,1,7,8,5)
      disease.list[[each.disease]] <- disease.list[[each.disease]][c(9,3,10,11,13,8,2,5,7,1,4,6,12)]
    } else {
      #order(all.frac.loh)
      #med 0.20 c(3,6,2,9,7,8,5,4,1)
      #med.seg.stroma c(10,2,7,8,1,4,6,9,5,3)
      #mu 0.25 c(9, 2, 1, 8, 7, 3, 6, 5, 4)
      disease.list[[each.disease]] <- disease.list[[each.disease]][c(10, 9, 2, 1, 8, 7, 3, 6, 5, 4)]
    }
  }
}

layout(matrix(c(rep(rep(0,30),padding.col),    # Padding whitespace on left for y-titles
                rep(c(0,45,45,45,seq(from=1, to=44, by=2),0,0,0,0), num.samples),    # Sample plots
                0,46,46,46,seq(from=2, to=44, by=2),0,0,0,0,    # Total count plots
                rep(rep(0,30),padding.col)  # Padding whitespace on right for y-titles  
                ), ncol=(2 + num.samples + 1 + 2), byrow=FALSE))




all.chr.int.df <- matrix(ncol=(4 + length(names(disease.list[[each.disease]]))), nrow=0)
colnames(all.chr.int.df) <- c("chrom", "Start.bp", "End.bp", names(disease.list[[each.disease]]), "t_loh")

for(each.chr in gsub("^", "chr", c(1:22))){
  print(paste("Dealing with ", each.disease, ":", each.chr, sep=""))
  
  ### Goes through each sample and creates an Interval list for each chromosome
  chr.matrix.list <- list() # Master Intervals list
  for(each.sample in names(disease.list[[each.disease]])){
    sample.df <- disease.list[[each.disease]][[each.sample]]
    sample.df <- as.data.frame(sample.df)
    
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
      sample.df <- as.data.frame(sample.df)
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
  
  if(mk.plots){
    num.of.dis <- length(names(chr.matrix.list))
    #  PLOT 1 - Plots the LOH fragment for EACH chromosome and EACH sample
    par(mar=c(0,4.1,0,0))
    plot(0, type='n', 
         xlim=c(0.5,(num.of.dis + 0.5)), ylim=c(1,ucsc.chr.row$size), 
         axes=FALSE, lty=0,
         ylab='', xlab='')
    #Create a function to generate a continuous color palette
    dis.cnt <- 1
    for(each.dis in names(chr.matrix.list)){  
      each.dis.df <- total.int.df[which(total.int.df[,each.dis] >= 1),]
      if(dim(each.dis.df)[1] > 0){
        dis.cn.col <- each.dis.df[,each.dis]
        if(any(dis.cn.col > 3)){
          dis.cn.col[which(dis.cn.col > 3)] <- 3
        }
        for(each.cn.t in c(0:3)){
          dis.cn.col <- gsub(paste("^", each.cn.t, "$", sep=""), dis.cn.pal[each.cn.t + 1], dis.cn.col) 
        }
        
        rect(xleft = (dis.cnt - 0.5), ybottom = each.dis.df$Start.bp,
             xright = (dis.cnt + 0.5), ytop = each.dis.df$End.bp,
             col=dis.cn.col, border=NA)  
        dis.cn.col <- NULL
      }
  
      dis.cnt <- dis.cnt + 1
    }
  
    # Separators for chromosomes
    #axis(1, at=c(0,(num.of.dis + 1)), labels=FALSE, tick=TRUE, col=alpha("grey", 0.50))
    #axis(3, at=c(0,(num.of.dis + 1)), labels=FALSE, tick=TRUE, col=alpha("grey", 0.50))
    #abline(v = seq(from=0.5, to=(dis.cnt - 0.5), by = 1))
    axis(side = 2, at=(ucsc.chr.row$size / 2), labels = gsub("^chr", "", each.chr), las=2, lty=0)
    if(each.chr == 'chr1'){
      axis(side = 2, at=c(ucsc.chr.row$size), labels = "", las=2) #Top  
    }
    axis(side = 2, at=c(0), labels = "", las=2) #Bottom
    if(each.chr == 'chr22'){
      axis(side = 1, at = c(1:(dis.cnt - 1)), labels = names(chr.matrix.list), las=2)
    }
    
    #  PLOT 2 - Plots the LOH fragment representing the aggregation of all samples
    par(mar=c(0,0,0,2.1))
    plot(0, type='n', 
         xlim=c(1,1), ylim=c(1,ucsc.chr.row$size),
         xaxt='n', yaxt='n', lty=0,
         ylab='', xlab='')
    #Create a function to generate a continuous color palette
    disPal <- colorRampPalette(c('white','gray22'))
    disPal <- disPal(num.of.dis + 1)
    dis.col <- total.int.df$t_loh
    for(each.loh.cnt in c(0:num.of.dis)){
      dis.col <- gsub(each.loh.cnt, disPal[each.loh.cnt + 1], dis.col) 
    }
    
    
    if(dim(total.int.df)[1] > 0){
      # TEMP Fix for not enough colours in a gradient
      if(length(dis.col > 7)) dis.col <- rep("white", length(total.int.df$Start.bp))
      rect(xleft = 0.5, ybottom = total.int.df$Start.bp,
           xright = 1.5, ytop = total.int.df$End.bp,
           col=dis.col, border=NA)
    }
    
    if(each.chr == 'chr22'){
      axis(side = 1, at = 1, labels = 'Total LOH', las=2)
    }
  }
    
  all.chr.int.df <- rbind(all.chr.int.df, total.int.df)
}

# Calculates the percent of the total genome that is LOH
all.chr.int.df[,'length'] <- all.chr.int.df[,'End.bp'] - all.chr.int.df[,'Start.bp']
all.frac.loh <- c()
for(each.sample.name in names(disease.list[[each.disease]])){
  print(each.sample.name)
  sample.df <- all.chr.int.df[which(all.chr.int.df[,each.sample.name] > 0),]
  sample.df <- as.data.frame(sample.df)
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

dev.off()
  

#save(disease.list, total.genome.size, ucsc.chrom, file="data.aggregateLoh.Rdata")
