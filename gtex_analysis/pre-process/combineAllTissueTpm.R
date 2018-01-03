args <- commandArgs(TRUE)
tpm.dir <- '/mnt/work1/users/pughlab/external_data/GTEx/FPKMmatrix/tpm/rdata'

all.tissue.tpm <- list.files(tpm.dir, pattern="tpm.Rdata")
all.tpm.mat <- data.frame()
for(each.tissue.tpm in all.tissue.tpm){
  load(file.path(tpm.dir, each.tissue.tpm))
  if(each.tissue.tpm %in% all.tissue.tpm[1]){
    all.tpm.mat <- tpm.mat
  } else if(rownames(tpm.mat) == rownames(all.tpm.mat)) {
    all.tpm.mat <- cbind(all.tpm.mat, tpm.mat)
  } else {
    print(paste("Warning: Could not append ", each.tissue.tpm, " - mismatching rows", sep=""))
  }
}

save(all.tpm.mat, file=file.path(tpm.dir, "allTissue.tpm.Rdata"))