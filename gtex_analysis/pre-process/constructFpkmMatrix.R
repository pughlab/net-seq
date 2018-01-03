args <- commandArgs(TRUE)
output.dir <- "/mnt/work1/users/pughlab/external_data/GTEx/FPKMmatrix"

tissue.dir <- args[1]
tissue.name <- args[2]
setwd(tissue.dir)


for(each.dir in list.files(tissue.dir, pattern="^G")){
  print(each.dir)
  if(file.exists(file.path(getwd(), each.dir, "genes.fpkm_tracking"))){
    
    sample.df <- read.table(file.path(getwd(), each.dir, "genes.fpkm_tracking"),
                            sep="\t", header=TRUE, stringsAsFactors=FALSE,
                            check.names=FALSE)
    sample.df <- sample.df[with(sample.df, order(tracking_id)),]
    
    # Parse the dataframe for FPKM
    sample.anno.df <- sample.df[,c("tracking_id", "gene_id", "gene_short_name", "locus")]
    sample.fpkm.df <- sample.df[,c("FPKM")]
    
    #Grab the annotation and dataframe structure from the first fpkm file
    if(each.dir %in% list.files(tissue.dir, pattern="^G")[1]){
      global.anno.df <<- sample.anno.df
      all.fpkm.df <<- data.frame(matrix(nrow=dim(global.anno.df)[1], ncol=0))
      rownames(all.fpkm.df) <- global.anno.df$tracking_id
      all.fpkm.df <- cbind(all.fpkm.df, sample.fpkm.df)
      colnames(all.fpkm.df)[dim(all.fpkm.df)[2]] <- each.dir
      # Else, check if the annotation is the same and append if it is
    } else {
      if(dim(sample.anno.df)[1] != dim(global.anno.df)[1]){
        print(paste("Warning: Dataframes are NOT the same size", each.dir, sep=""))
        print("Warning: Trying to appoint mismatched IDs...")
        all.fpkm.df$temp <- sample.df[match(rownames(all.fpkm.df), sample.df[,"tracking_id"]),"FPKM"]
        colnames(all.fpkm.df)[dim(all.fpkm.df)[2]] <- each.dir
      } else if (all(sample.anno.df == global.anno.df)){
        print(paste("Appended: ", each.dir, sep=""))
        all.fpkm.df <- cbind(all.fpkm.df, sample.fpkm.df)
        colnames(all.fpkm.df)[dim(all.fpkm.df)[2]] <- each.dir
      } else {
        print(paste("Warning: Trying to appoint mismatched IDs for ", each.dir, sep=""))
        all.fpkm.df$temp <- sample.df[match(rownames(all.fpkm.df), sample.df[,"tracking_id"]),"FPKM"]
        colnames(all.fpkm.df)[dim(all.fpkm.df)[2]] <- each.dir
      }
    }
  } else {
    print(paste("Ignored: No fpkm file found for: ", each.dir, sep=""))
  }
  
}

print(paste("     Dimension of fpkm dataframe: ", dim(all.fpkm.df)[2], sep=""))
print(paste("     Number of files processed: ", length(list.files(tissue.dir, pattern="^G")), sep=""))
print(paste("Saving to ", output.dir, sep=""))

save(all.fpkm.df, file=file.path(output.dir, paste(tissue.name, "fpkm.Rdata", sep=".")))
save(global.anno.df, file=file.path(output.dir, paste(tissue.name, "annotation.Rdata", sep=".")))

write.table(all.fpkm.df, file=file.path(output.dir, paste(tissue.name, "fpkm.tsv", sep=".")),
            quote=FALSE, sep="\t", row.names=TRUE, col.names=TRUE)
write.table(global.anno.df, file=file.path(output.dir, paste(tissue.name, "annotation.tsv", sep=".")),
            quote=FALSE, sep="\t", row.names=FALSE, col.names=TRUE)