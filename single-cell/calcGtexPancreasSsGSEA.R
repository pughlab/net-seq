###########################################################
## Create GeneSetCollection Object for Pancreas Genesets
## Date: November, 16, 2018
## By: Rene Quevedo
###########################################################

## Creates the Gene Set for GSEA for the specific cell types of Pancreas
library ("GSEABase")
library ("GSVA")
library ("GSVAdata")
library("org.Hs.eg.db")
data.dir <- file.path("git", "net-seq", "gtex_analysis", "data")
pancrea.cells <- read.table(file.path(data.dir, '20genes-cells.csv'), 
                            header=TRUE, stringsAsFactors = FALSE, sep=",")
ens.pancreas <- apply(pancrea.cells, 2, annotate)
pancreas.gs <- lapply(colnames(ens.pancreas), function(cell.type){
  GeneSet(setName=cell.type, 
          geneIdType = ENSEMBLIdentifier(),
          geneIds=unique(as.character(ens.pancreas[,cell.type, drop=TRUE])),
          collectionType = ExpressionSetCollection(),
          shortDescription = "20 Gene set that are descriptive of each cell type from the single cell pancreas atlas",
          longDescription = "",
          organism = "Homo Sapiens",
          pubMedIds = "27693023",
          urls = "http://www.ncbi.nlm.nih.gov/pubmed/27693023",
          contributor = "Muraro M.J.")
})
Pancreas.gsc <- GeneSetCollection(pancreas.gs)
save(Pancreas.gsc, 
     file = file.path("git", "net-seq", "gtex_analysis", "data", "ssPancreas.Rdata"))




#### Applies the Pancreas single-cell type to the GTEx dataset
#### LIBRARIES ####
library ("GSEABase")
library ("GSVA")
library ("GSVAdata")

annotate <- function(hgnc){
  require("biomaRt")
  require("org.Hs.eg.db")
  
  ens <- mapIds(org.Hs.eg.db,
                keys=hgnc,
                column="ENSEMBL",
                keytype="SYMBOL",
                multiVals="first")
  
  if(any(is.na(ens))){
    mart <- useMart("ENSEMBL_MART_ENSEMBL")
    mart <- useDataset("hsapiens_gene_ensembl", mart)
    na.ens <- getBM(mart=mart, 
                    attributes=c("external_gene_name", "hgnc_symbol", "ensembl_gene_id"),
                    filter="external_gene_name",
                    values=names(ens[which(is.na(ens))]),
                    uniqueRows=TRUE)
    ids <- sapply(na.ens$external_gene_name, function(x) {
      grep(paste0("^", x, "$"), names(ens), ignore.case=TRUE)
    })
    ens[unlist(ids)] <- rep(na.ens$ensembl_gene_id,
                            sapply(ids, length))
  }
  
  ens
}

gtex.dir <- "~/pughlab/external_data/GTEx/FPKMmatrix/tpm/tsv"
gtex.tpm.files <- list.files(gtex.dir, pattern = "Pancrea.*tpm.tsv$")
expr.data <- read.table(paste(gtex.dir,gtex.tpm.files[1],sep="/"), sep = "\t", header = T, stringsAsFactors = F);
row.names(expr.data) <- sapply(strsplit(row.names(expr.data), split = ".", fixed = T), function(x) return(x[1]));
# Log2 transform TPM values and change to matrix 
expr.data <- as.matrix(log2(expr.data+1));

load("/mnt/work1/users/pughlab/projects/NET-SEQ/rna_seq_external/clustering/ref/ssPancreas.Rdata")
expr.tiss <- expr.data[,1:5]
pancreas.scores <- gsva(expr.data[,1:5],
                        Pancreas.gsc,
                        method = "gsva",
                        min.sz = 2,max.sz = 500,
                        ssgsea.norm = T)








