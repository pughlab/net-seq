PDIR='/mnt/work1/users/pughlab/projects/NET-SEQ/rna_seq_external/MAE/samples'
setwd(PDIR)

genes <- c("MEN1", "ATRX", "DAXX")
MAD.samples <- lapply(paste0("samples_", genes, ".txt"), function(i) {
  samples <- read.table(i, header=FALSE,
                        stringsAsFactors = FALSE)
  samples$x <- 1
  samples 
})
names(MAD.samples) <- genes
metadata.raw <- read.table(file.path("attributes", "SraRunTable.txt"),
                       header = TRUE, stringsAsFactors = FALSE, 
                       sep="\t", fill=NA)
metadata <- metadata.raw[,c("Run", "pannet_subtype")]
colnames(metadata) <- c("TRACK_ID", "pannet_subtype")

MAD <- Reduce(function(x,y) merge(x, y, by="V1", all=TRUE), MAD.samples)
MAD[is.na(MAD)] <- 0
colnames(MAD) <- c("TRACK_ID", genes)
MAD$TRACK_ID <- gsub("\\..*", "", MAD$TRACK_ID)

MAD <- merge(MAD, metadata, by="TRACK_ID", all=TRUE)
write.table(MAD, file=file.path("attributes", "attributes.MAD.txt"),
            col.names = TRUE, row.names = FALSE,
            quote=FALSE, sep="\t")


save(MAD.samples, metadata.raw, MAD, file=file.path("attributes", "MAD.RData"))