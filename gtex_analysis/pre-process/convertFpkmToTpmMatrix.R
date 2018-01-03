fpkm.rdata.dir <- '/mnt/work1/users/pughlab/external_data/GTEx/FPKMmatrix/fpkm/rdata'
tpm.dir <- '/mnt/work1/users/pughlab/external_data/GTEx/FPKMmatrix/tpm'
rnaseq.tools <- '/mnt/work1/users/pughlab/src/R_wrappers/RNAseq_Processing/calculateTPMFromCounts.R'

dir.create(file.path(tpm.dir, "rdata"), recursive = TRUE)
dir.create(file.path(tpm.dir, "tsv"), recursive = TRUE)
source(rnaseq.tools)

all.tissues <- unique(gsub("\\..+", "", list.files(fpkm.rdata.dir, pattern="Rdata")))
for(each.tissue in all.tissues){
  load(file.path(fpkm.rdata.dir, paste(each.tissue, "fpkm.Rdata", sep=".")))
  load(file.path(fpkm.rdata.dir, paste(each.tissue, "annotation.Rdata", sep=".")))
  
  # convert FPKM to TPM
  print(paste("Generating TPM matrix for ", each.tissue, sep=""))
  tpm.mat <- calculateTPMFromFPKM(as.matrix(all.fpkm.df))
  
  if(all(dim(tpm.mat) == dim(all.fpkm.df)) &
       all(rownames(tpm.mat) == global.anno.df[,"tracking_id"])){
    save(tpm.mat, file=file.path(tpm.dir, "rdata", paste(each.tissue, "tpm.Rdata", sep=".")))
    save(global.anno.df, file=file.path(tpm.dir, "rdata", paste(each.tissue, "annotation.Rdata", sep=".")))
    
    write.table(tpm.mat, file=file.path(tpm.dir, "tsv", paste(each.tissue, "tpm.tsv", sep=".")),
                quote=FALSE, sep="\t", row.names=TRUE, col.names=TRUE)
    write.table(global.anno.df, file=file.path(tpm.dir, "tsv", paste(each.tissue, "annotation.tsv", sep=".")),
                quote=FALSE, sep="\t", row.names=TRUE, col.names=TRUE)
  } else {
    print(paste("Could not generate the TP matrix for: ", each.tissue, sep=""))
  } 
}
