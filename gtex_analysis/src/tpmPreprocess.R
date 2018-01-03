#GTEx TPM Functions

# Returns the tracking_id for the given gene_ID: Ensemble or HUGO
#   tracking_id: the rowname values for the TPM matrices
#   Uses the global.anno.df to convert a given gene.name (HUGO or Ensemble)
#     into the tracking_id for rowname retrieval
getId <- function(gene.name, id.type='hugo', global.anno.df=NA, ...){
  gene.tracking_id <- NULL
  if(id.type %in% 'hugo'){
    gene.row <- which(global.anno.df$gene_short_name %in% gene.name)
    gene.tracking_id <- global.anno.df[gene.row, 'tracking_id']
  } else if(id.type %in% 'ensemble'){
    #Convert ensemble gene ID from ENSG00000273367.10 -> ENSG00000273367
    ens.id <- gsub("\\..+$", "", global.anno.df$gene_id)
    ens.id <- gsub("\\..+$", "", gene.name)
    
    gene.row <- which(ens.id %in% gene.name)
    gene.tracking_id <- global.anno.df[gene.row, 'tracking_id']
  } else {
    print("Warning: unrecognized id.type.  Use 'hugo' or 'ensemble'")
  }
  
  return(gene.tracking_id)
}

# Gets the 95 CI for a vector of TPMs
get95Ci <- function(tpm.expr){
  mean.tpm <- mean(tpm.expr, na.rm=TRUE)
  sd.tpm <- sd(tpm.expr, na.rm=TRUE)
  
  error <- qt(0.975, df=(length(tpm.expr) - 1)) * sd.tpm/sqrt(length(tpm.expr))
  left <- mean.tpm - error
  right <- mean.tpm + error
  return(c(mean=mean.tpm, leftErr=left, rightErr=right))
  #return(paste(round(mean.tpm,1), "[", round(left, 1), ", ", round(right,1), "]", sep=""))
}

# Gets the tissue or subtype TPM value for a single gene from the TPM matrices created on Mordor
#  Returns a list where each element is the TPM matrix for that given type/subtype
# OR: if a gene name is provided, it will return a list of TPM values from each tissue/subtype for
#  that given gene
splitTissueTpm <- function(all.tissues, use.subtype=TRUE, gene.id=NA, id.type='hugo',
                           tpm.dir='/mnt/work1/users/pughlab/external_data/GTEx/FPKMmatrix/tpm/rdata', ...){
  # Format the annotations for easy extraction
  hist.types.list <- split(sra.dict, sra.dict$histology)
  
  # Get the TPM for the given gene name from each tissue
  tpm.tissue.list <- list()
  for(tissue.id in all.tissues){
    print(paste("Processing ", tissue.id, sep=""))
    load(file.path(tpm.dir, paste(tissue.id, ".annotation.Rdata", sep="")))
    load(file.path(tpm.dir, paste(tissue.id, ".tpm.Rdata", sep="")))
    
    # Retrieves tissue subtypes (e.g. Aorta) based on the annotation rather than the global 
    # tissue type (e.g. BloodVessel)
    if(use.subtype){
      # Get histology subtype specific TPMs (i.e. BloodVessel: Artery - Aorta, BloodVessel: Artery - Tibial)
      hist.subtypes.list <- split(hist.types.list[[tissue.id]], hist.types.list[[tissue.id]]$body.site)
      for(each.subtype in names(hist.subtypes.list)){
        print(paste("  subtype: ", each.subtype, sep=""))
        sample.index <- which(colnames(tpm.mat) %in% hist.subtypes.list[[each.subtype]]$bam.names)
        if(is.na(gene.id)){
          subtype.tpm <- tpm.mat[, sample.index]
        } else {
          gene.id <- getId(gene.name, id.type=id.type, global.anno.df=global.anno.df)
          subtype.tpm <- tpm.mat[gene.id, sample.index]
        }
        
        tpm.tissue.list[[each.subtype]] <- subtype.tpm
      }
    } else {
      if(is.na(gene.id)){
        tpm.tissue.list[[tissue.id]] <- tpm.mat[,]
      } else{
        tpm.tissue.list[[tissue.id]] <- tpm.mat[getId(gene.name, id.type=id.type, global.anno.df),]
      }
    }
  }
  return(tpm.tissue.list)
}


# Gets the correlation of a list of ordered TPM expressions (x) with all tissue sutbype samples (ordered)
getTissueCor <- function(x, tpm.list, cor.type){
  lapply(tpm.list, function(each.tissue.mat) {
    sapply(as.data.frame(each.tissue.mat), function(each.tissue.sample) {
      cor(x, each.tissue.sample, use="pairwise.complete.obs", method=cor.type) })
  })
}