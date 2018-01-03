####################################
###### clusterSampleAndPancan
######   Purpose: To cluster the PanCancer Segment file based on LOH
######    on a gene-wise basis.  This is done through identifying
######    LOH status for each gene and identifying % identity match
######    using a binary identifier.  Hierarchial clustering based on
######    percent LOH-identity match
####################################
library(intervals)
library(scales)
###############################
#         Variables
###############################
# Load pan-cancer ordered list derived from pancancer_loh.R
load("/Users/rquevedo/git/loh_cn_visualize/data/pancan_ordered_list.Rdata")

# PanCancer (ABSOLUTE) segment file
abs.seg.file <- '/Users/rquevedo/Desktop/PughLab/NET-seq/absolute_pan_can/seg_meta_data.csv'
ucsc.chrom.info.file <- '/Users/rquevedo/Desktop/bhk_lab/reference/ucsc.hg19.chromInfo.txt'     #hg19
remove.sex <- 1   # Removes the chr X and Y

#Sample Info
pnet.seg.file <- '/Users/rquevedo/Desktop/PughLab/NET-seq/Sequenza/bkup/all/segments/PNET_segments.txt'
sample.id <- 'pancreatic NET'
sinet.seg.file <- '/Users/rquevedo/Desktop/PughLab/NET-seq/Sequenza/bkup/all/segments/SINET_segments.txt'
sample.id2 <- 'small-intestine NET'
pnet.ccl.file <- '/Users/rquevedo/Desktop/PughLab/NET-seq/vandamme_cell_lines/adjusted_afplotter-loh_segments/PNET-CCL_segments.txt'
sample.id3 <- 'pNET cell lines'

#UCSC Refseq genes
refseq.f <- '/Users/rquevedo/Desktop/bhk_lab/reference/ucsc.refgene.txt'


###############################
#         Functions
###############################
getChrBp <- function(ucsc.chrom, x){
  x <- as.integer(as.character(x))
  chr.lengths <- ucsc.chrom[which(ucsc.chrom$chr.num < x), 'size']
  chr.size <- sum(as.numeric(as.character(chr.lengths)))
  return(chr.size)
}

addChrRect <- function(ucsc.chrom, num.of.samples, total.genome.size){
  mid.chr <- NULL
  for(each.row in 1:dim(ucsc.chrom)[1]){
    if(each.row != dim(ucsc.chrom)[1]){
      start.chr.bp <- ucsc.chrom[each.row,'chr.st.size']
      end.chr.bp <- ucsc.chrom[(each.row + 1),'chr.st.size']
    } else {
      start.chr.bp <- ucsc.chrom[each.row,'chr.st.size']
      end.chr.bp <- total.genome.size
    }
    
    rect(start.chr.bp, 0,
         end.chr.bp, (num.of.samples+ 1), 
         col=alpha("white", 0), border="black")
    mid <- (start.chr.bp + ((end.chr.bp - start.chr.bp) / 2))
    mid.chr <- c(mid.chr, mid)
  }
  return(mid.chr)
}

plotLohSegments <- function(x, col, ucsc.chrom){
  chr.num <- as.integer(as.character(x[['chr.num']]))
  sbp <- as.integer(as.character(x[['Start.bp']]))
  ebp <- as.integer(as.character(x[['End.bp']]))
  
  start.bp <- sbp + ucsc.chrom[match(chr.num, ucsc.chrom$chr.num),'chr.st.size']
  end.bp <- ebp + ucsc.chrom[match(chr.num, ucsc.chrom$chr.num),'chr.st.size']
  
  rect(start.bp, (sample.cnt - 0.25), 
       end.bp, (sample.cnt + 0.25), col=col, border=NA)
}

###############################
#           Main
###############################
# Absolute Segments
abs.seg <- read.csv(abs.seg.file,header = TRUE, sep = ",", quote = "\"", check.names=FALSE)
abs.seg <- abs.seg[which(abs.seg$TUMOR_TYPE != 'Cell-line'),]
abs.seg$Sample <- as.character(abs.seg$Sample)  #To deal with Factors for splitting
abs.seg$PRIMARY_DISEASE <- as.character(abs.seg$PRIMARY_DISEASE)  #To deal with Factors for splitting
abs.disease.seg.list <- split(abs.seg, abs.seg$PRIMARY_DISEASE)
abs.seg.list <- lapply(abs.disease.seg.list, function(x) split(x, x$Sample))

ordered.dis <- ordered.dis[! ordered.dis %in% c("pancreatic NET", "small-intestine NET")]
ordered.dis <- rev(ordered.dis)

#NET Segments
pnet.seg <- read.csv(pnet.seg.file, header = TRUE, sep = "\t", quote = "\"", check.names=FALSE)
pnet.seg.list <- split(pnet.seg, pnet.seg$Sample)
pnet.seg.list[['Sample']] <- NULL
sinet.seg <- read.csv(sinet.seg.file, header = TRUE, sep = "\t", quote = "\"", check.names=FALSE)
sinet.seg.list <- split(sinet.seg, sinet.seg$Sample)
sinet.seg.list[['Sample']] <- NULL
pnet.ccl.seg <- read.csv(pnet.ccl.file, header = TRUE, sep = "\t", quote = "\"", check.names=FALSE)
pnet.ccl.seg.list <- split(pnet.ccl.seg, pnet.ccl.seg$ID)
pnet.ccl.seg.list[['ID']] <- NULL

# Refseq Gene File
refgene.df <- read.csv(refseq.f, sep="\t", header=T)

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


###### Retrieves the LOH Genes:
# for(each.seg.list in list(sinet.seg.list, pnet.seg.list, abs.seg.list)){
#   for(each.seg in names(each.seg.list)){
#     print(paste("Retrieving genes for ", each.seg, "...", sep=""))
#     each.seg.list[[each.seg]] <- getSegmentGenes(each.seg.list[[each.seg]], refgene.df, pos.names=c(chr='Chromosome', start='Start.bp', end='End.bp'))
#   }
# }
# abs.seg.list2 <- lapply(abs.seg.list, function(x) getSegmentGenes(x, refgene.df, pos.names=c(chr='Chromosome', start='Start.bp', end='End.bp')))
# pnet.seg.list2 <- lapply(pnet.seg.list, function(x) getSegmentGenes(x, refgene.df, pos.names=c(chr='Chromosome', start='Start.bp', end='End.bp')))
# sinet.seg.list2 <- lapply(sinet.seg.list, function(x) getSegmentGenes(x, refgene.df, pos.names=c(chr='Chromosome', start='Start.bp', end='End.bp')))
pdf("loh_pancan_raw.pdf")
layout(matrix(c(0,2,2,2,2,
                0,2,2,2,2,
                0,2,2,2,2,
                0,1,1,1,1), nrow=4, byrow=TRUE))

par(mar=c(5.1,4.1,0.5,2.1))
names(pnet.seg.list) <- gsub("-N|absFormat/", "", names(pnet.seg.list))
names(sinet.seg.list) <- gsub("-N|absFormat/", "", names(sinet.seg.list))
if(length(grep("T_1", names(pnet.seg.list)))>0){
  pnet.seg.list <- pnet.seg.list[paste0("NET-00", c("1-T_2", "3-T_1", "3-T_2", "8-T_2", "9-T_1", "9-T_2"))]
  sinet.seg.list <- sinet.seg.list[paste0("NET-00", c("2-T_1", "2-T_2", "6-T_1", "6-T_2", "7-T_1", "7-T_2"))]
}
net.list <- list(pnets=pnet.seg.list, 
                 sinets=sinet.seg.list,
                 ccl=pnet.ccl.seg.list)
num.of.samples <- sum(unlist(lapply(net.list, function(x) length(names(x)))))
disease.list <- net.list
save(total.genome.size, disease.list, ucsc.chrom, file="~/git/loh_cn_visualize/data/data.aggregateLoh.T2-1.Rdata")

#generateLohPlot(net.list, num.of.samples, ucsc.chrom, total.genome.size)  # Function doesnt work, need to step through function at the moment

par(mar=c(0.5,4.1,4.1,2.1))
num.of.samples <- sum(unlist(lapply(abs.seg.list, function(x) length(names(x)))))
disease.list <- abs.seg.list
#generateLohPlot(abs.seg.list, num.of.samples, ucsc.chrom, total.genome.size) # Function doesnt work, need to step through function at the moment





generateLohPlot <- function(disease.list, num.of.samples, ucsc.chrom, total.genome.size){
  plot(0, type="n", xlim=c(1,total.genome.size), ylim=c(0,num.of.samples), 
       yaxt='n', xaxt='n', bty='n', xlab="", ylab="")  # Create a blank plot
  mid.chr <- addChrRect(ucsc.chrom, num.of.samples, total.genome.size)   # Add transparent boxes to show chromosomal edges
  axis(3, at=mid.chr, labels=ucsc.chrom$chrom, col="white", cex.axis=0.50, las=1)
  
  sample.cnt <- 1     # Used to indicate the boundaries of different diseases
  disease.col <- c("grey", "white")
  disease.cnt <- 1    # Used to color the different diseases
  disease.mids <- NULL
  #for(each.name in ordered.dis){                  # \   For PanCancer ordered list
  #  each.list <- disease.list[[each.name]]        # /
  for(each.list in disease.list){          # For disease samples
    pre.sample.cnt <- sample.cnt
    post.sample.cnt <- ( sample.cnt + length(names(each.list)) )
    
    # Add the Disease Sample alternating colours
    rect(1, (pre.sample.cnt - 0.35), 
         total.genome.size, (post.sample.cnt + 0.35), 
         col=alpha(disease.col[(disease.cnt %% 2) + 1], 0.5), 
         border=NA)
    d.mid <- pre.sample.cnt + ((post.sample.cnt - pre.sample.cnt) / 2)
    disease.mids <- c(disease.mids, d.mid)
  
    #for(each.sample in as.character(all.frac.loh.list[[each.name]]$Sample)){  # For Pancancer
    for(each.sample in names(each.list)){    # For our own disease list
      seg.df <- each.list[[each.sample]]
      seg.df.loh <- seg.df[which(seg.df$LOH == 1),]
      seg.df.loh$chr.num <- as.integer(as.character(gsub("^chr", "", seg.df.loh$Chromosome, ignore.case = TRUE)))
      if(dim(seg.df.loh)[1] > 0){
        for(each.loh in 1:dim(seg.df.loh)[1]){
          x <- seg.df.loh[each.loh,]  
        
          chr.num <- as.integer(as.character(x[['chr.num']]))
          sbp <- as.integer(as.character(x[['Start.bp']]))
          ebp <- as.integer(as.character(x[['End.bp']]))
          
          tcn.x <- as.integer(as.character(x$modal_A1)) + as.integer(as.character(x$modal_A2))

          if(tcn.x > 2){
            col="red"
          } else if(tcn.x < 2){
            col="dodgerblue"
          } else {
            col="thistle"
          }
          
          start.bp <- sbp + ucsc.chrom[match(chr.num, ucsc.chrom$chr.num),'chr.st.size']
          end.bp <- ebp + ucsc.chrom[match(chr.num, ucsc.chrom$chr.num),'chr.st.size']
          
          
          rect(start.bp, (sample.cnt - 0.25), 
               end.bp, (sample.cnt + 0.25), col=col, border=NA)
        }
      }
      sample.cnt <- sample.cnt + 1
    }
    
    disease.cnt <- disease.cnt + 1
  }
  axis(2, at=disease.mids, labels=names(disease.list), col="white", cex.axis=0.50, las=2)
  
}
dev.off()







