PDIR='/mnt/work1/users/pughlab/projects/NET-SEQ/rna_seq_external/gtex/sampleToGtexCor'
net <- file.path(PDIR, "input", "NETSEQ_on-treatment.tpm.pnet.Rdata")
srr <- file.path(PDIR, "input", "SRR.tpm.Rdata")
adm.ref <- file.path(PDIR, "ref", "attributes.MAD.txt")

load(net); net.tpm <- as.data.frame(tpm.mat)
load(srr); srr.tpm <- as.data.frame(tpm.mat)
adm <- read.table(adm.ref, sep="\t", header=TRUE, stringsAsFactors = FALSE)

net.tpm$ids <- rownames(net.tpm)
srr.tpm$ids <- rownames(srr.tpm)
pnet.tpm <- merge(net.tpm, srr.tpm, by="ids", all=TRUE)
rownames(pnet.tpm) <- pnet.tpm$ids
pnet.tpm <- pnet.tpm[,-grep("ids", colnames(pnet.tpm))]


mad.split <- split(adm, f=adm$pannet_subtype)
pnet.mad.tpm <- lapply(mad.split, function(x){
  idx <- match(x$TRACK_ID, colnames(pnet.tpm))
  as.matrix(pnet.tpm[,na.omit(idx)])
})

tpm.mat <- pnet.mad.tpm[[1]]
save(tpm.mat, file=file.path(PDIR, "input", "pnet-mutMAD.Rdata"))
tpm.mat <- pnet.mad.tpm[[2]]
save(tpm.mat, file=file.path(PDIR, "input", "pnet-wtMAD.Rdata"))



