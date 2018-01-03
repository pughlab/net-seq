.libPaths(c(.libPaths(), "/mnt/work1/users/pughlab/bin/r_lib/library/3.1"))
###############################################################
#
#                 Pan-Cancer LOH Plots
#
# Author: Rene Quevedo
# Date Created: Sep-16-2015
###############################################################
# Function: Creates a plot that compares all pan-cancer TCGA
#   samples analyzed by Absolute done by Carter at al. (2012) 
#   and plots the genomic fragment of LOH, CN-LOH, and aneuploidy
#   LOH in comparison to a given sample (NET)
###############################################################
library(scales) #for alpha colouring

###############################
#         Variables
###############################
# Old original dir/samples:
pnet.seg.file <- '~/git/net-seq/pancan_loh/data/PNET_segments.txt'
sample.id <- 'PNET (Discovery)'
sinet.seg.file <- '~/git/net-seq/pancan_loh/data/SINET_segments.txt'
sample.id2 <- 'GINET (Discovery)'
abs.seg.file <- '~/git/net-seq/pancan_loh/data/seg_meta_data.csv'
ucsc.chrom.info.file <- '~/git/net-seq/pancan_loh/data/ucsc.hg19.chromInfo.txt'     #hg19


#Pan-Cancer
remove.sex <- 1   # Removes the chr X and Y
sample.cutoff <- 10
alpha.val <- 0.6
plot.size <- 50
rm.disease.groups <- c("Pediatric GIST", "MFH", "Rhabdoid")


#Shallow WGS
pnet.otb.seg.file <- '~/git/net-seq/pancan_loh/data/pnet.wgsShallowCnv.sequenza.txt'
sample.otb.id <- 'PNET (Validation)'
sinet.otb.seg.file <- '~/git/net-seq/pancan_loh/data/ginet.wgsShallowCnv.sequenza.txt'
sample.otb.id2 <- 'GINET (Validation)'
#Samples to remove upon review from our resident neuroendocrine pathologist Dr. Sylvia Asa
pathology.blacklist <- paste("NET-2-0", c("01a", "01b", "08", "12", "14", "26", "28", "36"), sep="")  

###############################
#         Functions
###############################
{
  # Function: getLohInformation
  # Purpose:  Given a seg from abs, returns all LOH information to plot
  # Input:  x <- a row in abs.seg containing the headers "Sample, LOH, modal_A1, modal_A2, length, TUMOR_TYPE, PLATFORM"
  # Returns:  List containing the following information:
  #             tumor_type : Vector(char) : type of tumor the sample is from
  #             platform_type : Vector(char) : Array platform
  #             sample_name : Vector(char) : name of the sample
  #             cn_loh_seg : data.frame : contains all copy-neutral loh segments
  #             del_loh_seg : data.frame : contains all loh segments with TCN = 1
  #             aneu_loh_seg : data.frame : contains all loh segments with tCN >= 3
  #             tot_loh_seg : data.frame : contains all LOH segments
  #             loh_frac : data.frame : contains all genomic fraction and fraction of all LOH
  getLohInformation <- function(sample.seg){
    cn.loh.rows <- which(sample.seg$LOH == 1 & (as.numeric(as.character(sample.seg$modal_A1)) + as.numeric(as.character(sample.seg$modal_A2)) == 2))
    del.loh.rows <- which(sample.seg$LOH == 1 & (as.numeric(as.character(sample.seg$modal_A1)) + as.numeric(as.character(sample.seg$modal_A2)) == 1))
    aneu.loh.rows <- which(sample.seg$LOH == 1 & (as.numeric(as.character(sample.seg$modal_A1)) + as.numeric(as.character(sample.seg$modal_A2)) >= 3))
    tot.loh.rows <- which(sample.seg$LOH == 1)
    
    loh.frac.df <- data.frame(cn_loh = sum(as.numeric(as.character(sample.seg[cn.loh.rows,'length']))) / total.genome.size,
                              del_loh = sum(as.numeric(as.character(sample.seg[del.loh.rows,'length']))) / total.genome.size,
                              aneu_loh = sum(as.numeric(as.character(sample.seg[aneu.loh.rows,'length']))) / total.genome.size,
                              tot_loh = sum(as.numeric(as.character(sample.seg[tot.loh.rows,'length']))) / total.genome.size)
    loh.frac.df <- rbind(loh.frac.df, (loh.frac.df[1,] / loh.frac.df[1,'tot_loh']))
    rownames(loh.frac.df) <- c("frac_of_genome", "frac_of_tot_loh")
    
    
    sample.seg.list <- list(tumor_type = unique(sample.seg$TUMOR_TYPE),
                            platform_type = unique(sample.seg$PLATFORM),
                            primary_disease = unique(sample.seg$PRIMARY_DISEASE),
                            sample_name = unique(sample.seg$Sample),
                            cn_loh_seg = sample.seg[cn.loh.rows,],
                            del_loh_seg = sample.seg[del.loh.rows,],
                            aneu_loh_seg = sample.seg[aneu.loh.rows,],
                            tot_loh_seg = sample.seg[tot.loh.rows,],
                            loh_frac = loh.frac.df)
    return(sample.seg.list)
    
  }
  
  # Same as getLohInformation, but strictly looks for any deviation from modal_A1 and modal_A2 being 1 and 1 each (diploid)
  getCnInformation <- function(sample.seg){
    tot.cn.rows <- which(sample.seg$modal_A1 != 1 | sample.seg$modal_A2 != 1)
    gain.cn.rows <- which(as.numeric(as.character(sample.seg$modal_A1)) + as.numeric(as.character(sample.seg$modal_A2)) > 2)
    del.cn.rows <- which(as.numeric(as.character(sample.seg$modal_A1)) + as.numeric(as.character(sample.seg$modal_A2)) < 2)
    cn.loh.rows <- which(sample.seg$LOH == 1 & (as.numeric(as.character(sample.seg$modal_A1)) + as.numeric(as.character(sample.seg$modal_A2)) == 2))
    
    
    cn.frac.df <- data.frame(cn_loh = sum(as.numeric(as.character(sample.seg[cn.loh.rows,'length']))) / total.genome.size,
                              del_cn = sum(as.numeric(as.character(sample.seg[del.cn.rows,'length']))) / total.genome.size,
                              gain_cn = sum(as.numeric(as.character(sample.seg[gain.cn.rows,'length']))) / total.genome.size,
                              tot_cn = sum(as.numeric(as.character(sample.seg[tot.cn.rows,'length']))) / total.genome.size)
    cn.frac.df <- rbind(cn.frac.df, (cn.frac.df[1,] / cn.frac.df[1,'tot_cn']))
    rownames(cn.frac.df) <- c("frac_of_genome", "frac_of_tot_cn")
    
    
    sample.seg.list <- list(tumor_type = unique(sample.seg$TUMOR_TYPE),
                            platform_type = unique(sample.seg$PLATFORM),
                            primary_disease = unique(sample.seg$PRIMARY_DISEASE),
                            sample_name = unique(sample.seg$Sample),
                            cn_loh_seg = sample.seg[cn.loh.rows,],
                            del_cn_seg = sample.seg[del.cn.rows,],
                            gain_cn_seg = sample.seg[gain.cn.rows,],
                            tot_cn_seg = sample.seg[tot.cn.rows,],
                            cn_frac = cn.frac.df)
    return(sample.seg.list)
    
  }
  
  getLohRatio <- function(loh, tot_loh){
    loh <- as.numeric(as.character(loh))
    tot_loh <- as.numeric(as.character(tot_loh))
    loh.ratio <- loh/tot_loh
    loh.ratio[is.nan(loh.ratio)] <- 0
    
    mean(loh.ratio)
  }
  
  # Get X-coordinates for ordered loh fractions
  getDisXCoord <- function(x, plot.size, disease.order, boxplot.col, ord.col){
    x.start <- ((match(x$disease[1], disease.order) - 1) * plot.size) + 1
    x.end <- x.start + plot.size
    
    # Order based on tot_loh, cn_loh, aneu_loh, del_loh
    x <- x[order(x[,ord.col]),]
    
    # If there's only one sample
    if(dim(x)[1]==1) {
      x$X.pos<-x.start+(plot.size/2)
    } else {
      x$X.pos<-seq(from = x.start, 
                   to = x.end, 
                   by = plot.size / dim(x)[1])[-1]
    }
    
    #Add the colour for each disease
    x$col <- rep(boxplot.col[match(x$disease[1], disease.order)], dim(x)[1])
    
    return(x)
  }
  
  # Adds the Median Line to the LOH Plots
  plotMedLine <- function(x, row.num, mid.plot.range, plot.size, boxplot.col, lwd=2){
    spacing <- (plot.size / 2) - ((plot.size * 0.40) / 2)
    x0 <- mid.plot.range[row.num] - spacing
    x1 <- mid.plot.range[row.num] + spacing 
    y <- x[3]
    segments(x0, y, x1, y, col="black", lwd=2)
  }
  
  # Adds the confidence intervals to the LOH plots
  plotCiIntervals <- function(x, row.num, mid.plot.range, plot.size){
    spacing <- (plot.size / 2) - ((plot.size * 0.60) / 2)
    x0 <- mid.plot.range[row.num] - spacing
    x1 <- mid.plot.range[row.num] + spacing 
    segments(x0, x['low.ci'], x1, x['low.ci'], col=alpha("black", 0.50))
    segments(x0, x['high.ci'], x1, x['high.ci'], col=alpha("black", 0.50))
    segments(mid.plot.range[row.num], x['low.ci'], mid.plot.range[row.num], x['high.ci'], col=alpha("black", 0.50))
  }
  
  # Plots the LOH Points
  plotLohPoints <- function(loh.df, col.name, plot.col, alpha.val=0.5){
    loh.df[grep("grey", loh.df$col),'col'] <- plot.col
    
    points(loh.df$X.pos,
           as.numeric(as.character(loh.df[,col.name])),
           col=alpha(loh.df$col, alpha.val),
           pch=20, cex=0.5)
  }
  
  # Returns the confidence interval for a univariate t-distribution
  getConfInt <- function(x, col.name, conf.val){
    conf.val <- conf.val + (( 1 - conf.val) / 2 ) # e.g. 0.95
    col.num <- which(colnames(x) %in% col.name)  # e.g. "tot_loh"
    
    # Error value: (t-value | confidence) * (standard error of mean)
    err.val <- qt(conf.val, df=(dim(x)[1] - 1)) * (sd(x[,col.num])/sqrt(dim(x)[1]))
    # Calculate X confidence limits
    lower.ci <- mean(x[,col.num]) - err.val
    higher.ci <- mean(x[,col.num]) + err.val
    
    return(c(low.ci=lower.ci, high.ci=higher.ci))
  }
}
###############################
#           Main
###############################

# Read in absolute segment files formatted to include metadata
abs.seg <- read.csv(abs.seg.file,header = TRUE, sep = ",", quote = "\"", check.names=FALSE)
abs.seg <- abs.seg[which(abs.seg$TUMOR_TYPE != 'Cell-line'),]
abs.seg$Sample <- as.character(abs.seg$Sample)
abs.seg.list <- split(abs.seg, abs.seg$Sample)
#Discovery cohort
pnet.seg <- read.csv(pnet.seg.file, header = TRUE, sep = "\t", quote = "\"", check.names=FALSE)
pnet.seg.list <- split(pnet.seg, pnet.seg$Sample)
sinet.seg <- read.csv(sinet.seg.file, header = TRUE, sep = "\t", quote = "\"", check.names=FALSE)
sinet.seg.list <- split(sinet.seg, sinet.seg$Sample)
# Validation cohort
pnet.otb.seg <- read.csv(pnet.otb.seg.file, header = TRUE, sep = "\t", quote = "\"", check.names=FALSE)
pnet.otb.seg.list <- split(pnet.otb.seg, pnet.otb.seg$Sample)
sinet.otb.seg <- read.csv(sinet.otb.seg.file, header = TRUE, sep = "\t", quote = "\"", check.names=FALSE)
sinet.otb.seg.list <- split(sinet.otb.seg, sinet.otb.seg$Sample)



# read and filter ucsc.chrom.info file
ucsc.chrom <- read.csv(ucsc.chrom.info.file, header=TRUE, sep="\t", check.names=FALSE)
colnames(ucsc.chrom) <- gsub("#", "", colnames(ucsc.chrom))
ucsc.chrom <-  ucsc.chrom[regexpr("_", ucsc.chrom$chrom) == -1,]
ucsc.chrom <- ucsc.chrom[regexpr("^chrm$", ucsc.chrom$chrom, ignore.case = TRUE, perl=TRUE) == -1,]
if(remove.sex == 1){
  ucsc.chrom <- ucsc.chrom[regexpr("^chr[xy]$", ucsc.chrom$chrom, ignore.case = TRUE, perl=TRUE) == -1,]
}
total.genome.size <- sum(as.numeric(ucsc.chrom$size))


##### calculate the genomic fraction of LOH and TCN for pan-cancer samples
sample.seg.list <- lapply(abs.seg.list, function(x) getLohInformation(x))
all.diseases <- unlist(lapply(sample.seg.list, function(x) x[['primary_disease']]))
disease.loh <- list()
all.frac.loh.list <- list()
for(each.disease in unique(all.diseases)){
  genome.frac <- list(cn_loh = sapply(sample.seg.list[which(all.diseases == each.disease)], function(x) x[['loh_frac']]['frac_of_genome', 'cn_loh']),
                     del_loh = sapply(sample.seg.list[which(all.diseases == each.disease)], function(x) x[['loh_frac']]['frac_of_genome', 'del_loh']),
                     aneu_loh = sapply(sample.seg.list[which(all.diseases == each.disease)], function(x) x[['loh_frac']]['frac_of_genome', 'aneu_loh']),
                     tot_loh = sapply(sample.seg.list[which(all.diseases == each.disease)], function(x) x[['loh_frac']]['frac_of_genome', 'tot_loh']))
  disease.frac.loh.df <- do.call("cbind", genome.frac)
  all.frac.loh.list[[as.character(each.disease)]] <- as.data.frame(cbind(disease.frac.loh.df, disease=rep(each.disease, dim(disease.frac.loh.df)[1])))
  
  tot.loh.frac <- list(cn_loh = sapply(sample.seg.list[which(all.diseases == each.disease)], function(x) x[['loh_frac']]['frac_of_tot_loh', 'cn_loh']),
                       del_loh = sapply(sample.seg.list[which(all.diseases == each.disease)], function(x) x[['loh_frac']]['frac_of_tot_loh', 'del_loh']),
                       aneu_loh = sapply(sample.seg.list[which(all.diseases == each.disease)], function(x) x[['loh_frac']]['frac_of_tot_loh', 'aneu_loh']))
  disease.loh[[as.character(each.disease)]] <- list(genome_fraction = genome.frac,
                                                    total_loh_fraction = tot.loh.frac)
}

##### Calculate Genomic fraction of LOH and TCN for NET samples
#NET pre-processing for absolute calls
sinet.seg.list[['Sample']] <- NULL
pnet.seg.list[['Sample']] <- NULL
sinet.otb.seg.list[['Sample']] <- NULL
pnet.otb.seg.list[['Sample']] <- NULL
sinet.otb.seg.list[pathology.blacklist] <- NULL
pnet.otb.seg.list[pathology.blacklist] <- NULL

pnet.list <- lapply(pnet.seg.list, function(x) getLohInformation(x))
pnet.genome.frac <- list(cn_loh = sapply(pnet.list, function(x) x[['loh_frac']]['frac_of_genome', 'cn_loh']),
                        del_loh = sapply(pnet.list, function(x) x[['loh_frac']]['frac_of_genome', 'del_loh']),
                        aneu_loh = sapply(pnet.list, function(x) x[['loh_frac']]['frac_of_genome', 'aneu_loh']),
                        tot_loh = sapply(pnet.list, function(x) x[['loh_frac']]['frac_of_genome', 'tot_loh']))
pnet.frac.loh.df <- do.call("cbind", pnet.genome.frac)
pnet.frac.loh.df <- pnet.frac.loh.df[!rownames(pnet.frac.loh.df) %in% 'Sample',]

sinet.list <- lapply(sinet.seg.list, function(x) getLohInformation(x))
sinet.genome.frac <- list(cn_loh = sapply(sinet.list, function(x) x[['loh_frac']]['frac_of_genome', 'cn_loh']),
                         del_loh = sapply(sinet.list, function(x) x[['loh_frac']]['frac_of_genome', 'del_loh']),
                         aneu_loh = sapply(sinet.list, function(x) x[['loh_frac']]['frac_of_genome', 'aneu_loh']),
                         tot_loh = sapply(sinet.list, function(x) x[['loh_frac']]['frac_of_genome', 'tot_loh']))
sinet.frac.loh.df <- do.call("cbind", sinet.genome.frac)
sinet.frac.loh.df <- sinet.frac.loh.df[!rownames(sinet.frac.loh.df) %in% 'Sample',]

pnet.otb.list <- lapply(pnet.otb.seg.list, function(x) getLohInformation(x))
pnet.otb.genome.frac <- list(cn_loh = sapply(pnet.otb.list, function(x) x[['loh_frac']]['frac_of_genome', 'cn_loh']),
                         del_loh = sapply(pnet.otb.list, function(x) x[['loh_frac']]['frac_of_genome', 'del_loh']),
                         aneu_loh = sapply(pnet.otb.list, function(x) x[['loh_frac']]['frac_of_genome', 'aneu_loh']),
                         tot_loh = sapply(pnet.otb.list, function(x) x[['loh_frac']]['frac_of_genome', 'tot_loh']))
pnet.otb.frac.loh.df <- do.call("cbind", pnet.otb.genome.frac)
pnet.otb.frac.loh.df <- pnet.otb.frac.loh.df[!rownames(pnet.otb.frac.loh.df) %in% 'Sample',]

sinet.otb.list <- lapply(sinet.otb.seg.list, function(x) getLohInformation(x))
sinet.otb.genome.frac <- list(cn_loh = sapply(sinet.otb.list, function(x) x[['loh_frac']]['frac_of_genome', 'cn_loh']),
                          del_loh = sapply(sinet.otb.list, function(x) x[['loh_frac']]['frac_of_genome', 'del_loh']),
                          aneu_loh = sapply(sinet.otb.list, function(x) x[['loh_frac']]['frac_of_genome', 'aneu_loh']),
                          tot_loh = sapply(sinet.otb.list, function(x) x[['loh_frac']]['frac_of_genome', 'tot_loh']))
sinet.otb.frac.loh.df <- do.call("cbind", sinet.otb.genome.frac)
sinet.otb.frac.loh.df <- sinet.otb.frac.loh.df[!rownames(sinet.otb.frac.loh.df) %in% 'Sample',]


# -------- To filter for the 10 samples or more, re run from this poitn on and uncomment
# for(each.s in box.x.val$names[which(sample.size < sample.cutoff)]){
# all.frac.loh.list[[each.s]] <- NULL
# }

all.frac.loh.list[[sample.id]] <- as.data.frame(cbind(pnet.frac.loh.df, disease=rep(sample.id, dim(pnet.frac.loh.df)[1])))
all.frac.loh.list[[sample.id2]] <- as.data.frame(cbind(sinet.frac.loh.df, disease=rep(sample.id2, dim(sinet.frac.loh.df)[1])))
all.frac.loh.list[[sample.otb.id]] <- as.data.frame(cbind(pnet.otb.frac.loh.df, disease=rep(sample.otb.id, dim(pnet.otb.frac.loh.df)[1])))
all.frac.loh.list[[sample.otb.id2]] <- as.data.frame(cbind(sinet.otb.frac.loh.df, disease=rep(sample.otb.id2, dim(sinet.otb.frac.loh.df)[1])))



# Removal of trouble datasets (Pediatric GIST and MFH)
if(length(rm.disease.groups) >= 1){
  rm.elements <- which(names(all.frac.loh.list) %in% rm.disease.groups)
  all.frac.loh.list <- all.frac.loh.list[-rm.elements]
}

# Calculate the median for total_loh per disease subtype
disease.median <- sapply(all.frac.loh.list, function(x) median(as.numeric(as.character(x$tot_loh))))
all.frac.loh.df <- do.call("rbind", all.frac.loh.list[names(disease.median[order(disease.median)])])
all.frac.loh.df$tot_loh <- as.numeric(as.character(all.frac.loh.df$tot_loh))




# Create the fraction of LOH  to total_loh [OUTDATED TO CREATE STACKED BARPLOT]
disease.cn.mean <- sapply(split(all.frac.loh.df,  all.frac.loh.df$disease), function(x) getLohRatio(x$cn_loh, x$tot_loh))
disease.del.mean <- sapply(split(all.frac.loh.df,  all.frac.loh.df$disease), function(x) getLohRatio(x$del_loh, x$tot_loh))
disease.aneu.mean <- sapply(split(all.frac.loh.df,  all.frac.loh.df$disease), function(x) getLohRatio(x$aneu_loh, x$tot_loh))

t.tot.loh.frac <- t(data.frame(cn=disease.cn.mean,
                    del=disease.del.mean,
                    aneu=disease.aneu.mean))


###############################
#           Pre-Visualization

# Obtain calculations for the paper:
# x.df <- all.frac.loh.df[-which(all.frac.loh.df$disease %in% c(sample.id, sample.id2, sample.otb.id, sample.otb.id2)),]
# apply(x.df[,c(1,2,3,4)], 2, function(x) median(as.numeric(as.character(x)), na.rm=TRUE))
# t.test(x=all.frac.loh.df[which(all.frac.loh.df$disease %in% sample.id), 'tot_loh'],
#        y=all.frac.loh.df[which(all.frac.loh.df$disease %in% 'Lung squamous'), 'tot_loh'],  alternative = 'greater')
{
  # Sample sizes and identify Tumour Sample of Interest
  box.x.val <- boxplot(tot_loh~disease, data=all.frac.loh.df, ylim=c(0,1))
  sample.size <- box.x.val$n
  # sample.pos <- which(box.x.val$names %in% c(sample.id, sample.id2))    # DISCOVERY Cohort
  # sample.otb.pos <- which(box.x.val$names %in% c(sample.otb.id, sample.otb.id2))  #VALIDATION Cohort
  sample.pos <- which(box.x.val$names %in% c(sample.id, sample.otb.id))  #PNETS
  sample.otb.pos <- which(box.x.val$names %in% c(sample.id2, sample.otb.id2))  #GINETS
  
  
  # Create Colour Scheme
  boxplot.col <- rep(alpha("black", alpha.val), length(sample.size))
  boxplot.col[which(sample.size < sample.cutoff)] <- alpha("black", alpha.val)
  boxplot.col[sample.pos] <- 'red'
  boxplot.col[sample.otb.pos] <- 'blue'
  
  # Identify the X.pos and Colours for each sample
  all.frac.loh.df$X.pos <- rep(0, dim(all.frac.loh.df)[1])
  disease.order <- unique(all.frac.loh.df$disease)
  all.frac.loh.df <- do.call("rbind", lapply(split(all.frac.loh.df, all.frac.loh.df$disease), 
                                             function(x) getDisXCoord(x, plot.size, disease.order, boxplot.col, 'tot_loh')))
  aneu.frac.loh.df <- do.call("rbind", lapply(split(all.frac.loh.df, all.frac.loh.df$disease), 
                                             function(x) getDisXCoord(x, plot.size, disease.order, boxplot.col, 'aneu_loh')))
  aneu.frac.loh.df$aneu_loh <- as.numeric(as.character(aneu.frac.loh.df$aneu_loh))
  cn.frac.loh.df <- do.call("rbind", lapply(split(all.frac.loh.df, all.frac.loh.df$disease), 
                                             function(x) getDisXCoord(x, plot.size, disease.order, boxplot.col, 'cn_loh')))
  cn.frac.loh.df$cn_loh <- as.numeric(as.character(cn.frac.loh.df$cn_loh))
  del.frac.loh.df <- do.call("rbind", lapply(split(all.frac.loh.df, all.frac.loh.df$disease), 
                                             function(x) getDisXCoord(x, plot.size, disease.order, boxplot.col, 'del_loh')))
  del.frac.loh.df$del_loh <- as.numeric(as.character(del.frac.loh.df$del_loh))
  
  # Retrieve the 95% Confidence limits for total_loh
  conf.int.list <- lapply(split(all.frac.loh.df, all.frac.loh.df$disease), function(x) getConfInt(x, "tot_loh", 0.95))
}


#           Visualization
###############################

{
  dir.create("~/Desktop/netseq/genie/plots/", recursive = TRUE, showWarnings = FALSE)
  pdf("~/Desktop/netseq/genie/plots/loh_pancan_net.na.genfrac.prim.pdf")
  #layout(matrix(c(0,0,0,0,1,1,2,2,3,3,4,4,0), nrow=13, byrow=TRUE))
  layout.matrix <- matrix(c(0, 1, 0, 0.01,
                            0, 1, 0.01, 0.011,
                            0, 1, 0.011, 0.012,
                            0, 1, 0.012, 0.35, 
                            0, 1, 0.35, 0.7), byrow=TRUE, ncol=4)
  split.screen(layout.matrix)
  # PLOT1: Incremental Scatter Median plot for LOH Genomic Fraction
  screen(5)
  par(mar=c(0, 4.1, 0, 4.1))
  plot(0, type='n', xlim=c(1,(plot.size * length(disease.order))), ylim=c(0,1), xaxt='n', xlab="", 
       ylab="% genome LOH", yaxt="n")
  plotLohPoints(all.frac.loh.df, "tot_loh", "grey", alpha.val)
  
  # Plot median lines for each tumor type
  mid.plot.range <- seq(plot.size, plot.size * length(disease.order), by=plot.size) - (plot.size / 2)
  sapply(1:ncol(box.x.val$stats), function(x) plotMedLine(box.x.val$stats[,x], x, mid.plot.range, plot.size, boxplot.col))
  #sapply(1:length(conf.int.list), function(x) plotCiIntervals(conf.int.list[[x]], x, mid.plot.range, plot.size))
  
  axis(2, at=seq(0,1,by=0.2), labels=seq(0, 100, by=20), las=2)
  # Label the tumor name and sample size
  disease.names <- apply(data.frame(dis=disease.order, size=sample.size), 1, function(x) paste(x[1], " (", as.integer(x[2]), ")", sep=""))
  dis.names <- disease.names 
  dis.names[which(sample.size < sample.cutoff)] <- ''
  dis.names <- gsub("Acute lymphoblastic leukemia", "ALL", dis.names)
  axis(3, at=mid.plot.range, labels=dis.names, las=2, col.axis='black', cex.axis=0.7)
  dis.names <- disease.names
  dis.names[which(sample.size >= sample.cutoff)] <- ''
  axis(3, at=mid.plot.range, labels=dis.names, las=2, col.axis=alpha('black', (alpha.val+0.2)), cex.axis=0.7)
  dis.names <- disease.names
  dis.names[-sample.pos] <- ''
  axis(3, at=mid.plot.range, labels=dis.names, las=2, col.axis='red', cex.axis=0.7)
  dis.names <- disease.names
  dis.names[-sample.otb.pos] <- ''
  axis(3, at=mid.plot.range, labels=dis.names, las=2, col.axis='blue', cex.axis=0.7)
  
  
  # # PLOT 2/3/4: Scatterplot for LOH Segment Gains, Neutral, Losses
  # #GAINS
  # screen(3)
  # par(mar=c(0, 4.1, 0.5, 4.1))
  # plot(0, type='n', xlim=c(1,(plot.size * length(disease.order))), ylim=c(0,1), xaxt='n', xlab="", ylab="Gain", yaxt="n")
  # plotLohPoints(aneu.frac.loh.df, "aneu_loh", "red")
  # sapply(1:ncol(box.x.val$stats), function(x) plotMedLine(boxplot(aneu_loh~disease, data=aneu.frac.loh.df, plot=FALSE)$stats[,x],
  #                                                                 x, mid.plot.range, plot.size, boxplot.col))
  # axis(2, yaxp=c(0,1,2), las=2)
  
  #NEUTRAL
  screen(4)
  par(mar=c(0, 4.1, 0.5, 4.1))
  plot(0, type='n', xlim=c(1,(plot.size * length(disease.order))), ylim=c(0,1), xaxt='n', xlab="", 
       ylab="% genome copy-neutral", yaxt="n")
  plotLohPoints(cn.frac.loh.df, "cn_loh", "thistle", 0.8)
  sapply(1:ncol(box.x.val$stats), function(x) plotMedLine(boxplot(cn_loh~disease, data=cn.frac.loh.df, plot=FALSE)$stats[,x],
                                                                x, mid.plot.range, plot.size, boxplot.col))
  axis(2, at=seq(0,1,by=0.2), labels=seq(0, 100, by=20), las=2)
  
  # #LOSS
  # screen(2)
  # par(mar=c(0, 4.1, 0.5, 4.1))
  # plot(0, type='n', xlim=c(1,(plot.size * length(disease.order))), ylim=c(0,1), xaxt='n', xlab="", ylab="Loss", yaxt="n")
  # plotLohPoints(del.frac.loh.df, "del_loh", "dodgerblue")
  # sapply(1:ncol(box.x.val$stats), function(x) plotMedLine(boxplot(del_loh~disease, data=del.frac.loh.df, plot=FALSE)$stats[,x],
  #                                                                 x, mid.plot.range, plot.size, boxplot.col))
  # axis(2, yaxp=c(0,1,2), las=2)
  close.screen(all.screens=TRUE)
  dev.off()
}
