# Function: formatClinvarDf
# Purpose:  Subsets for the assembly and regex extraction of Ref and Alt allele for SNPs
# Input:  clinvar.df : A dataframe containing the clinvar csv download
#         assembly: GRCh37 or GRCh38
# Returns:  dataframe : essentiall clinvar.df with ref_ and alt_  allele subsetted for assembly
formatClinvarDf <- function(clinvar.df, assembly){
  #Subset for only the proper assembly
  clinvar.df <- clinvar.df[which(clinvar.df$Assembly == assembly),]
  
  # Extract the Allele information from the Name column (NM_014630.2(ZNF592):c.3136G>A (p.Gly1046Arg))
  clinvar.snps.m <- regexpr(pattern="[ACGT]+>[ACGT]+", text=clinvar.df$Name, perl=TRUE)
  clinvar.snps <- regmatches(m = clinvar.snps.m, x = clinvar.df$Name)
  clinvar.indels <- which(clinvar.snps.m == -1)
  #Substring for ref and alt allele
  clinvar.snp.div <- regexpr(pattern=">", clinvar.snps, fixed = TRUE)
  ref_allele <- substring(clinvar.snps, 1, (clinvar.snp.div - 1))
  alt_allele <- substring(clinvar.snps, (clinvar.snp.div + 1))
  
  clinvar.df[which(clinvar.snps.m != -1),'ref_allele'] <- ref_allele
  clinvar.df[which(clinvar.snps.m != -1),'alt_allele'] <- alt_allele
  
  return(clinvar.df)
}


# Function: generateGenomeBins
# Purpose:  Uses the UCSC ChromInfo file to generate bins of a specific size for later on
#   CN evaluation
# Input:  - chr.df <- ucsc.chrom.info dataframe containing Chr and Size
#         - bin.s <- The size of the bins you want to create
# Returns:  genome.binned: data-frame containing the genome binned into the set sizes
generateGenomeBins <- function(chr.df, bin.s){
  # Formats the ucsc chrom info dataframe in an easier to use way.
  chr.df <- chr.df[ chr.df$chrom %in% global.chrom,]
  rownames(chr.df) <- chr.df$chrom
  chr.df <- chr.df[,-1]
  
  # Creates a List containing binned intervals for ezch chromsome
  chr.bin.list <- list()
  for(each.chr in global.chrom){
    bin.start.pos <- 1
    bin.max.matrix <- matrix(ncol=2, nrow=0)
    while(bin.start.pos < max(chr.df[as.character(each.chr),'size'])){
      bin.end.pos <- bin.start.pos + bin.s
      if(bin.end.pos > max(chr.df[as.character(each.chr),'size'])){
        bin.end.pos <- max(chr.df[as.character(each.chr),'size'])
      }
      bin.max.matrix <- rbind(bin.max.matrix, c(bin.start.pos, bin.end.pos))
      bin.start.pos <- bin.end.pos + 1
    }
    chr.bin.list[[as.character(each.chr)]] <- bin.max.matrix
  }
  
  return(chr.bin.list)
}



# Function: generateWeightedCompositeBin
# Purpose:  Looks at predetermined bins and proportionally weighs each CN-segment that makes up that bin
#   to return a weighted-CN-segment for that bin.
# Input:  composite.bin.df <- dataframe containing the actual intervals that make up each predetermined bins
#         seg.col.name <- column name for the segment header (i.e. seg_mean, hscr.a1, hscr.a2, etc...)
# Returns:  composite.bin.final.df: data-frame containing the genome binned into the set sizes and having a 
#   proportionally weighted CN segment according to to the segment length in that bin.
generateWeightedCompositeBin <- function(composite.bin.df, seg.col.name){
  composite.bin.final.df <- as.data.frame(matrix(ncol=3, nrow=0))
  for(each.uniq.bin in unique(composite.bin.df$Start.bin)){
    unique.bin.df <- composite.bin.df[composite.bin.df$Start.bin == each.uniq.bin, ]
    
    start.bin <- each.uniq.bin
    end.bin <- unique.bin.df[1,'End.bin']
    
    # Calculates the proportionally-weighted CN segment for the given bins    
    unique.bin.df$length <- unique.bin.df$End.bp - unique.bin.df$Start.bp
    total.cov <- sum(unique.bin.df$length)
    unique.bin.df$weighted.seg <- (unique.bin.df$length / total.cov) * unique.bin.df[,seg.col.name]
    unique.bin.df$weighted.copy.ratio <- (unique.bin.df$length / total.cov) * unique.bin.df[,'copy.ratio']
    
    composite.sum <- sum(unique.bin.df$weighted.seg)
    composite.copy.ratio.sum <- sum(unique.bin.df$weighted.copy.ratio)
    composite.bin.final.df <- rbind(composite.bin.final.df, 
                                    data.frame(start.bin, end.bin, composite.sum, composite.copy.ratio.sum))
  }
  colnames(composite.bin.final.df) <- c("Start.bin", "End.bin", seg.col.name, 'copy.ratio')
  
  return(composite.bin.final.df)
}


# Function: populateGenomeBins
# Purpose:  Takes in a pre-determined blank list of all bin sizes for each chromosome in the human genome
#   and maps out the given segment length and CN-segments to each of the bins.  This will allow for an easier
#   bin-vs-bin comparison in future steps.
# Input:  cn.val.f1: dataframe containing the CN segments and their seg-values (hscr.a1, seg_mean, etc..)
#         chr.bins: The blank list containing all preformatted empty bins
#         seg.col.name: column name for the segment header (i.e. seg_mean, hscr.a1, hscr.a2, etc...)
# Returns:  genome.tbins.df: Dataframe containing all the preformatted bins with their associated seg-values 
#   according to cn.val.f1.
populateGenomeBins <- function(cn.val.f1, chr.bins, seg.col.name){
  genome.tbins.df <- data.frame(matrix(ncol=5, nrow=0))
  for(each.chr in global.chrom){
    composite.bin.df <- as.data.frame(matrix(ncol=6, nrow=0)) #Keeps track of segments that do not fit whole bins
    matching.intervals.df <- as.data.frame(matrix(ncol=4, nrow=0)) # Keeps track of whole bin segments.
    
    colnames(matching.intervals.df) <- c("Start.bin", "End.bin", seg.col.name, 'copy.ratio')
    blank.chr.int <- Intervals(chr.bins[[each.chr]], closed=TRUE)
    cn.chr.f1 <- cn.val.f1[cn.val.f1$Chromosome %in% each.chr,]
    for(each.chr.row in 1:dim(cn.chr.f1)[1]){
      cn.int.f1 <- Intervals(matrix(c(cn.chr.f1[each.chr.row, 'Start.bp'], cn.chr.f1[each.chr.row, 'End.bp']), ncol=2))
      match.matrix <- as.matrix(interval_intersection(cn.int.f1, blank.chr.int))
      # Matches the blank matrix to the "match matrix" that has removed first and last row (to be calculated later)
      temp.matrix <- matrix(sort(intersect(chr.bins[[each.chr]], match.matrix[-c(1,nrow(match.matrix)),])), ncol=2, byrow=T)
      temp.matrix <- cbind(temp.matrix, rep(cn.chr.f1[each.chr.row, seg.col.name], nrow(temp.matrix)))
      temp.matrix <- cbind(temp.matrix, rep(cn.chr.f1[each.chr.row, 'copy.ratio'], nrow(temp.matrix)))
      colnames(temp.matrix) <- c("Start.bin", "End.bin", seg.col.name, 'copy.ratio')
      matching.intervals.df <- rbind(matching.intervals.df, temp.matrix)
      colnames(matching.intervals.df) <- c("Start.bin", "End.bin", seg.col.name, 'copy.ratio')
      
      # Stores the intervals that do not fill a complete bin in the composite dataframe to handle after.
      # bp-level intervals match in the bin
      first.match.row <- match.matrix[1,]
      last.match.row <- match.matrix[nrow(match.matrix),]
      # the bin being matched to
      first.match.blank.int <- Intervals(first.match.row)
      last.match.blank.int <- Intervals(last.match.row)
      composite.matrix <- matrix(c(first.match.row, blank.chr.int[!is.na(interval_overlap(blank.chr.int, first.match.blank.int) > 0)],
                                   last.match.row, blank.chr.int[!is.na(interval_overlap(blank.chr.int, last.match.blank.int) > 0)]), 
                                 nrow=2, byrow=T)
      composite.matrix <- cbind(composite.matrix, 
                                rep(cn.chr.f1[each.chr.row, seg.col.name], nrow(composite.matrix)),
                                rep(cn.chr.f1[each.chr.row, 'copy.ratio'], nrow(composite.matrix)))
      if(length(which(first.match.row == last.match.row)) ==  2){   # If the first.match.row vector matches last.match.row vector
        composite.bin.df <- data.frame(rbind(as.matrix(composite.bin.df), composite.matrix[1,]))
      } else {
        composite.bin.df <- data.frame(rbind(as.matrix(composite.bin.df), composite.matrix))
      }
      # Composite matrix result: start.bp, end.bp, start.bin, end.bin, segment.col.name
      
    }
    colnames(composite.bin.df) <- c("Start.bp", "End.bp", "Start.bin", "End.bin", seg.col.name, 'copy.ratio')
    matching.intervals.df <- rbind(matching.intervals.df, generateWeightedCompositeBin(composite.bin.df, seg.col.name))
    matching.intervals.df <- matching.intervals.df[order(matching.intervals.df$Start.bin),]
    #Finds the intervals not included in the analysis and fills them in
    match.intervals <- Intervals(matrix(c(matching.intervals.df[, 'Start.bin'], 
                                          matching.intervals.df[, 'End.bin']), ncol=2, byrow=F))
    mismatched.int.matrix <- as.matrix(interval_difference(blank.chr.int, match.intervals))
    mismatched.int.matrix<- cbind(mismatched.int.matrix, 
                                  rep(mismatch.val, dim(mismatched.int.matrix)[1]), 
                                  rep(mismatch.val, dim(mismatched.int.matrix)[1]))
    matching.intervals.df <- data.frame(rbind(as.matrix(matching.intervals.df), mismatched.int.matrix))
    
    matching.intervals.df <- matching.intervals.df[order(matching.intervals.df$Start.bin),]
    matching.intervals.df$chr <- each.chr
    
    
    
    genome.tbins.df <- rbind (genome.tbins.df, matching.intervals.df)
    # Need to find the complement of the intervals using interval_complement
  }
  colnames(genome.tbins.df) <- c("Start.bin", "End.bin", seg.col.name, 'copy.ratio', "Chr")
  return(genome.tbins.df)
}
