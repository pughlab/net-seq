.libPaths(c(.libPaths(), "/mnt/work1/users/pughlab/bin/r_lib/library/3.1"))
library(intervals);
library(AnnotationDbi)
library(org.Hs.eg.db)
#########################
#   Variables
#########################
#Debug
refgene.file <- '/Users/rquevedo/git/cn-annotater/refData/refgene.exon.bed'
get.strand <- 0

# seg.file <- '/Users/rquevedo/git/cn-annotater/example_data/aacr_pnet_genie.seg';
# seg.df <- read.csv(file=seg.file, header=T, sep="\t",
#                    stringsAsFactors=FALSE, check.names=FALSE)

#########################
#   Functions
#########################




#colnames(seg.df) <- c("ID", "chrom", "loc.start", "loc.end", "num.mark", "seg.mean")

# Chromosome list

annotateSeg <- function(seg.df, chrom.ids=NA, 
                        loc.headers=c("chrom", "loc.start", "loc.end")){
  print("Reading in the refgene file...")
  refgene.df <- read.csv(file=refgene.file, header=T, sep="\t", check.names=FALSE)
  seg.df <- validateColnames(seg.df, loc.headers)

  match.df <- data.frame(matrix(nrow=0, ncol=4))
  colnames(match.df) <- c('chrom', 'loc.start', 'loc.end', 'gene')
  if(is.na(chrom.ids)) global.chrom.ids <- as.vector(unique(seg.df$chrom))
  print(global.chrom.ids)
  
  for(chr.id in global.chrom.ids){
    print(paste("Retrieving genes for ", chr.id, "...", sep=""))
    # Adjust for the difference between chr1 and 1 for example
    if(length(grep("^chr", chr.id)) == 0){
      refgene.chr.id <- paste("chr", chr.id, sep="")
    } else {
      refgene.chr.id <- chr.id
    }
    
    #Subsets each dataset per chromosome and retrieves the location + CN value
    seg.chr.df <- subset(seg.df, chrom==chr.id, select=c("loc.start", "loc.end"))
    refgene.chr.df <- subset(refgene.df, chrom==refgene.chr.id, select=c("cdsStart", "cdsEnd", 'strand', 'name2'))
    
    #Creates an interval object based on the subsetted chromosomal segment start and ends
    seg.chr.interval <- Intervals(c(seg.chr.df$loc.start, seg.chr.df$loc.end))
    refgene.chr.interval <- Intervals(c(refgene.chr.df$cdsStart, refgene.chr.df$cdsEnd))
    
    # Overlap the segment intervals with the refgene intervals
    overlap.rows <- interval_overlap(seg.chr.interval, refgene.chr.interval)  
    overlap.rows.na <- lapply(overlap.rows, function(x) if(length(x) == 0) NA else x) #replace all integer(0) with NA to maintain order
    overlap.rows.gene <- lapply(overlap.rows.na, function(x) paste(unique(refgene.chr.df[x,'name2']), collapse = ","))
    seg.chr.df$gene <- unlist(overlap.rows.gene)
    
    
    seg.chr.df <- cbind(rep(chr.id, dim(seg.chr.df)[1]), seg.chr.df)
    colnames(seg.chr.df) <- c('chr', 'start', 'end', 'gene')
    match.df <- rbind(match.df, seg.chr.df)
  }
  
  if(dim(seg.df)[1] == dim(match.df)[1]){
    seg.df$gene <- match.df$gene
  }
  return(seg.df)
}


validateColnames <- function(seg.df, loc.headers){
  # Assumes loc.headers is composed of a 3-value vector
  #   [1] - chromosome identifier (e.g. chr)
  #   [2] - start position identifier (e.g. start.pos)
  #   [3] - end position identifier (e.g. end.pos)
  
  if(!all(loc.headers == c("chrom", "loc.start", "loc.end"))){
    chr.idx <- grep(loc.headers[1], colnames(seg.df))
    start.idx <- grep(loc.headers[2], colnames(seg.df))
    end.idx <- grep(loc.headers[3], colnames(seg.df))
    
    colnames(seg.df)[c(chr.idx, start.idx, end.idx)] <- c("chrom", "loc.start", "loc.end")
  }
  return(seg.df)
}
# 

# Attribute a copy-ratio/seg.mean to each gene, as well as annotate the Ensembl ID
mapGenes <- function(seg.df, col.id=NA){
  if(is.na(col.id)){
    gene.list <- apply(seg.df, 1, function(x) unlist(strsplit(x['gene'], split = ",")))
    names(gene.list) <- with(seg.df, paste(chrom, paste(round((loc.start / 1000000),2),  
                                                        round((loc.end / 1000000),2), sep="-"), 
                                           sep=":") )
  } else {
    seg.list <- split(seg.df, seg.df[,col.id])
    gene.list <- lapply(seg.list, function(x) unlist(strsplit(x$gene, ",")))
  }
  
  
  ord.genes <- unlist(gene.list)
  genes.df <- data.frame("seg.mean"=names(ord.genes), "genes"=ord.genes)
  genes.df$ensembl <- mapIds(org.Hs.eg.db,
                             keys=as.character(genes.df$genes),
                             keytype="SYMBOL",
                             column="ENSEMBL",
                             multiVals="first")
  return(genes.df)
}
# write.table(seg.df, file=output.file, col.names=TRUE,
#             quote=F, sep="\t", eol="\n", row.names=FALSE)
# write.table(genes.df, file=output.genes, col.names=TRUE,
#             quote=F, sep="\t", eol="\n", row.names=FALSE)