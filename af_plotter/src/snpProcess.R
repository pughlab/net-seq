# Takes a tumor-only mutect VCF and returns AF stats using variantAnnotation package
getMutectAf <- function(vcf, col.sel=1){
  print("** Notice: MuTect mode is only set up for single sample run at the moment")
  geno.gt <- geno(vcf)$GT[,col.sel]   # Genotype call
  geno.dp <- geno(vcf)$DP[,col.sel]   # Total depth
  geno.ad <- geno(vcf)$AD[,col.sel]   # Depth of A and B allele
  
  #Calculates the Allelic Fraction
  allele.cnt.df <- do.call("rbind", geno.ad)
  af.df <- allele.cnt.df / geno.dp  #convert to allelic fraction
  af.df <- cbind(af.df, allele.cnt.df, 
                 geno.dp, gsub(":.+", "", rownames(af.df)))
  colnames(af.df) <- c("a", "b", "raw_a", "raw_b", "depth", "chr")
  
  return(af.df)
}

# Takes a single-sample HaplotypeCaller VCF and returns AF stats using variantAnnotation package
getHaplotypeCallerAf <- function(vcf){
  # Parse genotype of VCF
  geno.gt <- geno(vcf)$GT   # Genotype call
  geno.dp <- geno(vcf)$DP   # Total depth
  geno.ad <- geno(vcf)$AD   # Depth of A and B allele
  
  #NET-001 Fix (due to merging of two samples and mixed headers, removes Sample_D1)
  if(sample.match.m %in% "NET-001-T_2.vcf.gz"){
    geno.gt <- geno.gt[,-which(colnames(geno.gt) %in% "Sample_D1"),drop=FALSE]
    geno.dp <- geno.dp[,-which(colnames(geno.dp) %in% "Sample_D1"),drop=FALSE]
    geno.ad <- geno.ad[,-which(colnames(geno.ad) %in% "Sample_D1"),drop=FALSE]
  }
  # / END Fix
  
  #Checks for two calls, an A and B allele (No AAB or ABB)
  allele.cnt <- apply(geno.ad, 1, function(x) unlist(x))   #unlist A and B allele calls
  num.of.alleles <- lapply(allele.cnt, function(x) length(x))    # Check how many calls (A and B calls)
  not.2.alleles <- which(unlist(num.of.alleles) != 2)
  only2.alleles <- which(!seq(1:length(allele.cnt)) %in% not.2.alleles)
  allele.cnt.only2 <- allele.cnt[only2.alleles]     # Only have 2 
  allele.cnt.only2.df <- do.call("rbind", allele.cnt.only2)
  
  #Calculates the Allelic Fraction
  af.df <- allele.cnt.only2.df / geno.dp[only2.alleles]  #convert to allelic fraction
  af.df <- cbind(af.df, allele.cnt.only2.df[,1], allele.cnt.only2.df[,2])  #append Allelic Count
  af.df <- cbind(af.df, geno.dp[only2.alleles,1])   # append Allelic Depth
  af.df <- cbind(af.df, gsub(":.+", "", rownames(af.df)))  # append Chr
  colnames(af.df) <- c("a", "b", "raw_a", "raw_b", "depth", "chr")
  
  return(af.df)
}

# Removes factors and adds genomic positions to AF dataframe
formatAfDf <- function(af.df){
  # Double check the formats and factors to ensure numeric and not character
  for(each.col in c("a", "b", "raw_a", "raw_b")){
    af.df[,each.col] <- as.numeric(as.character(af.df[,each.col]))
  }
  af.df <- as.data.frame(af.df)
  
  # Add genomic position
  m <- regexpr(":.+?_",rownames(af.df))
  chr.pos <- gsub("[:_]", "", regmatches(rownames(af.df), m))
  af.df[,'pos'] <- chr.pos
  return(af.df)
}

# Function: generateChrList
# Purpose:  Generates a dataframe from VCF file containing the Allelic Depth, Total Depth, and Genotype...
#   Separated by chromosome.
# Input:  sample.match.m : A vector value containing the VCF file name
#         vcf.dir: the directory the vcf file can be found in
# Returns:  List: Contains dataframes for each chromosome and the allelic depth and ratio
generateChrList <- function(sample.match.m, vcf.dir, col.id=1,
                            snp.caller='haplotypecaller', clinvar.stat=TRUE){
  require(VariantAnnotation)
  vcf <- readVcf(file.path(vcf.dir, sample.match.m), "hg19")
  z.af.list <- list()
  
  #### Sets up the Allelic Fraction dataframe
  # IF there is only one sample in the VCF file (HaplotypeCaller)
  if(snp.caller %in% 'haplotypecaller'){
    af.df <- getHaplotypeCallerAf(vcf)
    af.df <- formatAfDf(af.df)
  } else if(snp.caller %in% 'mutect'){
    af.df <- getMutectAf(vcf, col.sel=col.id)
    af.df <- formatAfDf(af.df)
  } else if(snp.caller %in% 'varscan2'){
    geno.rd <- geno(vcf)$RD[,2]  #Referencec Allele
    geno.ad <- geno(vcf)$AD[,2]  # Alternate Allele
    geno.dp <- geno(vcf)$DP[,2]  # Raw read depth
    
    #Obtain position location via regexp
    m <- regexpr(":.+?_",names(geno.ad))
    chr.pos <- gsub("[:_]", "", regmatches(names(geno.ad), m))
    #Obtain chr location
    chr.num <- gsub(":.+", "", names(geno.dp))
    
    af.df <- data.frame(a=(geno.rd/geno.dp),
                          b=(geno.ad/geno.dp),
                          depth=geno.dp,
                          chr=chr.num,
                          pos=chr.pos)
    
    #Remove rows in the dataframe that contain "NA" values
    na.rows <- apply(is.na(af.df), 1, any)
    af.df <- af.df[!na.rows,]
    
  }
  
  #Add in allele information
  #ref allele
  ref.allele.m <- regexpr("_[ACGT]+", rownames(af.df), ignore.case = FALSE, perl = TRUE)
  af.df$ref_allele <- substring(regmatches(x = rownames(af.df), m=ref.allele.m), 2)
  #alt allele
  alt.allele.m <- regexpr("/[ACGT]+", rownames(af.df), ignore.case = FALSE, perl = TRUE)
  af.df$alt_allele <- substring(regmatches(x = rownames(af.df), m=alt.allele.m), 2)
  #Remove certain Chr
  blacklist.chr <- c('chrM', 'chrY')
  
  #### Identifies Clinvar Variants
  af.df$clinvar_state <- rep("Unknown", dim(af.df)[1])
  af.df$clinvar_gene <- rep(NA, dim(af.df)[1])
  af.df$clinvar_gene_start <- rep(NA, dim(af.df)[1])
  af.df$clinvar_gene_end <- rep(NA, dim(af.df)[1])
  
  if(clinvar.stat){
    snp.match <- which((as.numeric(gsub("^chr", "", af.df$chr)) %in% clinvar.df$Chromosome) &
                         (as.numeric(af.df$pos) %in% clinvar.df$Start))
    
    # Function: getRefgenePos
    # Purpose:  Returns the CDS start and stop for a gene given it's name from Refgene
    getRefgenePos <- function(gene.name){
      gene.pos <- refgene.df[which(refgene.df$name2 %in% gene.name)[1],c('cdsStart', 'cdsEnd')]
      return(gene.pos)
    }
    
    for(each.clinvar.var in snp.match){
      snp.chr <- gsub("^chr", "", af.df[each.clinvar.var, 'chr'])
      snp.pos <- af.df[each.clinvar.var, 'pos']
      snp.ref <- af.df[each.clinvar.var, 'ref_allele']
      snp.alt <- af.df[each.clinvar.var, 'alt_allele']
      
      #Find all clinvar matches
      clinvar.m <- which(clinvar.df$Chromosome %in% snp.chr & 
                           clinvar.df$Start %in% snp.pos)
      
      if(length(clinvar.m) > 0){
        #Cycles through every possible option
        for(each.var in c(1:dim(clinvar.df[clinvar.m,])[1])){
          #checks the ref and alt alleles
          if(clinvar.df[clinvar.m[each.var],'ref_allele'] %in% snp.ref &
               clinvar.df[clinvar.m[each.var], 'alt_allele'] %in% snp.alt){          
            
            gene.pos <- getRefgenePos(as.character(clinvar.df[clinvar.m[each.var], 'GeneSymbol']))
            
            af.df[each.clinvar.var,'clinvar_gene_start'] <-  gene.pos$cdsStart
            af.df[each.clinvar.var,'clinvar_gene_end'] <-  gene.pos$cdsEnd
            af.df[each.clinvar.var,'clinvar_gene'] <-  as.character(clinvar.df[clinvar.m[each.var], 'GeneSymbol'])
            af.df[each.clinvar.var,'clinvar_state'] <-  as.character(clinvar.df[clinvar.m[each.var], 'ClinicalSignificance'])
          }
        }
      }
    }
  }
  
  # Splits into a list based on chromosomes
  if(length(grep("^chr", af.df$chr)) == 0) { af.df$chr <- gsub("^", "chr", af.df$chr) }
  chrSplit.af.list <- split(af.df, af.df$chr)
  blacklist.idx <- which(names(chrSplit.af.list) %in% blacklist.chr)
  if(length(blacklist.idx) > 0){
    chrSplit.af.list <- chrSplit.af.list[-blacklist.idx]
  }

  return(chrSplit.af.list)
}


# Function: overlapAf
# Purpose:  Looks at overlaps between RNA and Exome calls
#           Excludes filtering out SNV and SNP calls that are found only in one or the other
# Input:  ret.type : dna or rna
#         dna.list : DNA SNP calls from generateChrList
#         rna.list : RNA SNP calls from generateChrList
# Returns:  List: Contains dataframes for each chromosome and the allelic depth and ratio - overlap filtered
overlapAf <- function(ret.type, dna.list, rna.list){
  ret.list <- list()
  for(each.chr in names(dna.list)){
    dna.chr.list <- dna.list[[each.chr]]
    rna.chr.list <- rna.list[[each.chr]]
    
    # Filter for overlaps
    if(ret.type == 'dna'){
      ret.list[[each.chr]] <- mergeAf(dna.chr.list, rna.chr.list)
    } else if(ret.type == 'rna'){
      ret.list[[each.chr]] <- mergeAf(rna.chr.list, dna.chr.list)
    }
  }
  
  return(ret.list)
}

# Function: mergeAf
# Purpose:  Takes a single chr-level dataframe (or two), and separates pathological SNPs and SNVs from the df
#   prior to manipulating the df with: 
#       merge.stat = 1 : overlap dna and rna calls
#       merge.stat = min.depth : a minimum depth to filter out certain SNP calls that will cause noise
#   Reattaches the SNPs and SNVs afterwards.
# Input:  genome.chr.df.a : dna or rna snp df for chromosomes
#         genome.chr.df.b : dna or rna snp df for chromosomes, or NULL if merge.stat != 1
#         merge.stat : refer to Purpose
# Returns:  Data frame containing the filtered/manipulated dataframe with SNPs and SNVs intact
mergeAf <- function(genome.chr.df.a, genome.chr.df.b, merge.stat=1){
  # Partition out the SNVs and SNPs
  snp.df <- genome.chr.df.a[grep("Pathogenic", genome.chr.df.a$clinvar_state),]
  snv.df <- genome.chr.df.a[which(!is.na(genome.chr.df.a$Variant_Classification)),]
  
  if(merge.stat == 1){
    genome.chr.df.a <- genome.chr.df.a[which(rownames(genome.chr.df.a) %in% rownames(genome.chr.df.b)),]  
  } else{
    genome.chr.df.a <- genome.chr.df.a[which(genome.chr.df.a$depth > merge.stat),]
  }
  
  
  # Merge and remove duplicates
  t.genome.chr.df.a <- rbind(genome.chr.df.a, snv.df, snp.df)
  t.genome.chr.df.a <- t.genome.chr.df.a[order(as.numeric(as.character(t.genome.chr.df.a$pos))),]
  if(length(which(duplicated(t.genome.chr.df.a$pos))) > 0){
    t.genome.chr.df.a <- t.genome.chr.df.a[-which(duplicated(t.genome.chr.df.a$pos)),]
  }
  
  return(t.genome.chr.df.a)
}


# Segments based on a target depth
getShallowSnps <- function(x, targ.depth){
  x[which(as.integer(as.character(x$depth)) >= targ.depth),]
}

# Function: binSnps
# Purpose:  Takes a genome coord bin list and processed SNPs list (separated by chr)
#   and assigns each SNP to a genome coord bin, calculates the mean +/- 1 s.d.
# Inputs:   af.list: generateChrList processed snp list separated by chr
#           chrom.bins.list: genome coord bins list separated by chr
# Returns:  Similar format to af.list, separated by chr
binSnps <- function(af.list, shallow.wgs=FALSE, rm.hom=FALSE, ...){
  af.list[['chrX']] <- NULL
  orig.af.name <- names(af.list)
  af.names <- gsub("^chr", "", names(af.list), ignore.case=TRUE)
  names(af.list) <- af.names
  bin.names <- gsub("^chr", "", names(chrom.bins.list), ignore.case=TRUE)
  
  # If shallow WGS (1x-2x), filters for only 2-read depth SNPs
  if(shallow.wgs){
    af.list <- lapply(af.list, getShallowSnps, targ.depth=2)
  }

  
  bin.af.list <- lapply(af.names, function(each.chr){
    #Removes homozygous SNPs if true
    if(rm.hom){
      x.df <- af.list[[as.character(each.chr)]][which(as.numeric(as.character(af.list[[as.character(each.chr)]]$a)) > 0 &
                                                        as.numeric(as.character(af.list[[as.character(each.chr)]]$a)) < 1),]
      af.list[[as.character(each.chr)]] <- x.df
    }
    
    # Convert to interval object and look for which SNPs are in which genome bins
    snp.int <- Intervals(matrix(as.integer(c(af.list[[as.character(each.chr)]][,'pos'],
                                             af.list[[as.character(each.chr)]][,'pos'])),
                                byrow=FALSE, ncol=2), closed=TRUE)
    bin.int <- Intervals(chrom.bins.list[[as.character(each.chr)]], closed=TRUE)
    snp.bin.overlap <- interval_overlap(snp.int, bin.int)
    
    # Get the mean and +/- 1 sd of the binned SNPs
    
    if(length(unique(unlist(snp.bin.overlap))) == 1){
      snp.bins <- list()
      snp.bins[[1]] <- as.vector(sapply(unique(unlist(snp.bin.overlap)), function(x) which(snp.bin.overlap == x)))
    } else {
      snp.bins <- sapply(unique(unlist(snp.bin.overlap)), function(x) which(snp.bin.overlap == x))
    }
    # Given a dataframe from generateChrList format, will return the mean
    # and +/- 1 s.d. for all SNPs in that dataframe
    getBinStats <- function(x, allele.col, shallow.stat=FALSE){
      bin.x <- as.numeric(as.character(af.list[[as.character(each.chr)]][x, allele.col]))
  
      if(shallow.stat){
#         het.cnt <- length(bin.x[which(bin.x == 0.5)])
#         hom.cnt <- length(bin.x[which(bin.x != 0.5)])
        hom.cnt <- length(bin.x[which(bin.x == 1 | bin.x == 0)])
        het.cnt <- length(bin.x[which(bin.x < 1 & bin.x > 0)])
        return(c(het=het.cnt, hom=hom.cnt))
      } else {
        mean.x <- mean(bin.x)
        lower.sd.x <- mean.x - sd(bin.x)
        upper.sd.x <- mean.x + sd(bin.x)
        return(c("mean"=mean.x, "low.sd" = lower.sd.x, "high.sd" = upper.sd.x))
      }
    }

    if(length(snp.bins) == 1){
      chr.stat.df <- cbind(rep(gsub("^", "chr", each.chr), length(snp.bins)),    #Chr name repeated
                           t(chrom.bins.list[[as.character(each.chr)]][unique(unlist(snp.bin.overlap)),]),  # Pos of bins with SNPs
                           sapply(snp.bins, length))
    } else {
      chr.stat.df <- cbind(rep(gsub("^", "chr", each.chr), length(snp.bins)),    #Chr name repeated
                           chrom.bins.list[[as.character(each.chr)]][unique(unlist(snp.bin.overlap)),],  # Pos of bins with SNPs
                           sapply(snp.bins, length))
    }
    
    
    if(shallow.wgs){
      cnt.stat <- t(sapply(snp.bins, getBinStats, allele.col='a', shallow.stat=TRUE))
      chr.stat.df <- cbind(chr.stat.df,
                           cnt.stat)
      colnames(chr.stat.df) <- c("chr", "start.pos", "end.pos", "n",
                                 "het", "hom")
      chr.stat.df <- as.data.frame(chr.stat.df)
      
      #changes characters to integers and numerics
      chr.stat.df$start.pos <- as.integer(as.character(chr.stat.df$start.pos))
      chr.stat.df$end.pos <- as.integer(as.character(chr.stat.df$end.pos))
      chr.stat.df$n <- as.integer(as.character(chr.stat.df$n))
      chr.stat.df$het <- as.integer(as.character(chr.stat.df$het))
      chr.stat.df$hom <- as.integer(as.character(chr.stat.df$hom))
    } else {
      a.stat <- t(sapply(snp.bins, getBinStats, allele.col='a'))
      b.stat <- t(sapply(snp.bins, getBinStats, allele.col='b'))
      
      # Assembles the dataframe
      colnames(a.stat) <- paste(colnames(a.stat), ".a", sep="")
      colnames(b.stat) <- paste(colnames(b.stat), ".b", sep="")
      chr.stat.df <- cbind(chr.stat.df, 
                           a.stat,   #a-allele statistics
                           b.stat)   #b-allele statistics
      
      colnames(chr.stat.df) <- c("chr", "start.pos", "end.pos", "n",
                                 "mean.a", "low.sd.a", "high.sd.a",
                                 "mean.b", "low.sd.b", "high.sd.b")
      chr.stat.df <- as.data.frame(chr.stat.df)
      
      #changes characters to integers and numerics
      chr.stat.df$start.pos <- as.integer(as.character(chr.stat.df$start.pos))
      chr.stat.df$end.pos <- as.integer(as.character(chr.stat.df$end.pos))
      chr.stat.df$n <- as.integer(as.character(chr.stat.df$n))
      chr.stat.df$mean.a <- as.numeric(as.character(chr.stat.df$mean.a))
      chr.stat.df$low.sd.a <- as.numeric(as.character(chr.stat.df$low.sd.a))
      chr.stat.df$high.sd.a <- as.numeric(as.character(chr.stat.df$high.sd.a))
      chr.stat.df$mean.b <- as.numeric(as.character(chr.stat.df$mean.b))
      chr.stat.df$low.sd.b <- as.numeric(as.character(chr.stat.df$low.sd.b))
      chr.stat.df$high.sd.b <- as.numeric(as.character(chr.stat.df$high.sd.b))
    }
    
    return(chr.stat.df)
  })
  
  names(bin.af.list) <- orig.af.name
  return(bin.af.list)
}


# Filters out regions where there are a lot of SNPs in x-sized bin
# Returns the same each.chr.df but without the problem bin SNPs
#     default: 5snps in 1kb
filtRegion <- function(x, bin=5000, quantcut=0.995, default.ratio=5){
  require(intervals)
  #calculate number of breaks needed to get x bp break-sizes
  pos.spread <- max(x$pos) - min(x$pos)
  break.size <- ceiling(pos.spread / bin)
  
  # Find bins that have large amount of SNPs in them
  h.conf <- hist(x$pos, breaks=break.size, plot=FALSE)
  qcutoff <- quantile(h.conf$counts, quantcut)
  if(qcutoff < (default.ratio * (bin/1000))){
    print(paste("Changing minimum cutoff:  Calculated cutoff is too low - ", qcutoff, sep=""))
    qcutoff <- (default.ratio * (bin/1000))
  } else {
    print(paste("Filtering at ", qcutoff, sep=""))
  }
  
  # Find the trouble intervals
  hi.bins <- which(h.conf$counts > qcutoff)
  hi.bin.mat <- matrix(c(h.conf$breaks[hi.bins], 
                         h.conf$breaks[hi.bins + 1]), byrow=FALSE, ncol=2)
  snp.mat <- matrix(c(x$pos, x$pos), byrow=FALSE, ncol=2)
  hi.int <- Intervals(hi.bin.mat, closed=TRUE)
  snp.int <- Intervals(snp.mat, closed=TRUE)
  
  # Finds the SNPs that don't overlap with that region
  hi.overlap <- interval_overlap(snp.int, hi.int)
  x.df <- x[which(is.na(hi.overlap > 0)),]
  return(x.df)
}


