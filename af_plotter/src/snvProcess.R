# Function: annotateSnvs
# Purpose:  Uses mutect-oncotator calls to annotate mutations of interest
# Input:  maf.file : Name of the mutect-oncotator MAF file
#         maf.dir : Absolute path the mutect-oncotator MAF file directory
#         genome.list : DNA or RNA SNP calls from generateChrList
annotateSnvs <- function(maf.file, maf.dir, genome.list){
  input.maf <- read.csv(file.path(maf.dir, maf.file), header=TRUE, sep="\t", comment.char="#")
  input.maf.list <- split(input.maf, input.maf$Chromosome)
  
  ret.list <- list()
  for(each.chr in names(genome.list)){
    if(each.chr %in% gsub("^", "chr", names(input.maf.list))){
      # Merge Variants and HUGO for SNVs
      chr.merge.df <- merge(x = genome.list[[each.chr]], 
                            y = input.maf.list[[gsub("^chr", "", each.chr)]][,c("Hugo_Symbol", "Start_position", "Variant_Classification")],
                            by.x = "pos", by.y = "Start_position",
                            all.x = TRUE)
      chr.merge.df <- chr.merge.df[order(as.numeric(as.character(chr.merge.df$pos))),]
      rownames(chr.merge.df) <- rownames(genome.list[[each.chr]])
      # Replace to include the annotated data.frame
      ret.list[[each.chr]] <- chr.merge.df
    } else {
      # apply blank columns if no SNV is found on that chromosome
      chr.blank.df <- genome.list[[each.chr]]
      chr.blank.df$Hugo_Symbol <- NA
      chr.blank.df$Variant_Classification <- NA
      rownames(chr.blank.df) <- rownames(genome.list[[each.chr]])
      ret.list[[each.chr]] <- chr.blank.df
    }
    
  }
  return(ret.list)
}

removeHomSnps <- function(ret.type, dna.list, rna.list){
  ret.list <- list()
  for(each.chr in names(dna.list)){
    dna.chr.list <- dna.list[[each.chr]]
    rna.chr.list <- rna.list[[each.chr]]
    
    hom.pos <- as.numeric(as.character(dna.chr.list[which(dna.chr.list$a %in% c(0,1) & dna.chr.list$b %in% c(0,1)), 'pos']))
    if(length(hom.pos) > 0){
      # Filter for overlaps
      if(ret.type == 'dna'){
        ret.list[[each.chr]] <- dna.chr.list[-which(dna.chr.list$pos %in% hom.pos),]
      } else if(ret.type == 'rna'){
        ret.list[[each.chr]] <- rna.chr.list[-which(rna.chr.list$pos %in% hom.pos),]
      }
    } else {
      if(ret.type == 'dna'){
        ret.list[[each.chr]] <- dna.chr.list
      } else if(ret.type == 'rna'){
        ret.list[[each.chr]] <- rna.chr.list
      }
    }
    
  }
  return(ret.list)
}
