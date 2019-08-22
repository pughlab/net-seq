library(metafor)
setwd('~/git/net-seq/cgh_analysis/cgh_plotter/')
load("~/git/net-seq/cgh_analysis/data/cgh_df.RData")

specify_decimal <- function(x, k) format(round(x, k), nsmall=k)
getTableIndices <- function(table.idx){
  # plot ylims
  tbl.ylim <- (sum(table.idx) + 16)
  
  # plot Row indices
  group1.idx <- c(3:(table.idx[1] + 2))
  group2.idx <- c((max(group1.idx) + 6):(max(group1.idx) + 5 + table.idx[2]))
  group3.idx <- c((max(group2.idx) + 5):(max(group2.idx) + 4 + table.idx[3]))
  rows.idx <- c(group1.idx, group2.idx, group3.idx)
  
  #subgroup indices
  subgroup.idx <- c("G3"=max(group3.idx + 1),
                    "G2"=max(group2.idx + 1),
                    "G1"=max(group1.idx + 1))
  
  #column heading indices
  column.idx <- c((subgroup.idx["G3"] + 1),
                  (subgroup.idx["G2"] + 1),
                  (subgroup.idx["G1"] + 1),
                  "Title"=(as.integer(subgroup.idx["G3"]) + 3))
  
  #sumamry polygon indices
  polygon.idx <- c("G3"=min(group3.idx) - 1.5,
                   "G2"=min(group2.idx) - 1.5,
                   "G1"=min(group1.idx) - 1.5)
  
  
  return(list("ylim"=tbl.ylim,
              "rows"=rows.idx,
              "subs"=subgroup.idx,
              "cols"=column.idx,
              "poly"=polygon.idx))
}

#####  Gets the contigency table for each feature being studied
#####  and removes single confounders such as unknowns or the single LNET across all studies
cgh.contigency.list <- list()
metadata.clust <- c("Study.id", "NET", "F.stat", "Met.stat")
for (each.feat in metadata.clust){
  if (each.feat %in% 'NET'){
    temp.cgh.df <- cgh.df[which(!cgh.df$NET %in% 'LNET'),] #Remove the single LNET
    
    # Generate Contigency table and save data
    temp.table.list <- lapply(split(temp.cgh.df, temp.cgh.df$Study.id), function(x){
                          t(table(factor(x$NET, levels=c("PNET", "GINET")), 
                                  factor(x$ci.stat, levels=c("High-CI", "Low-CI"))))
                        })
    cgh.contigency.list[[each.feat]] <- do.call("rbind", lapply(temp.table.list, as.vector))
    cgh.contigency.list[[each.feat]] <- cbind(cgh.contigency.list[[each.feat]],
                                              rep(each.feat, dim(cgh.contigency.list[[each.feat]])[1]),
                                              rownames(cgh.contigency.list[[each.feat]]))
  }
  
  if(each.feat %in% 'F.stat'){
    temp.cgh.df <- cgh.df[which(!cgh.df[,'F.stat'] %in% 'Unk'),] #Remove the Unknown functional status samples
    
    # Generate Contigency table and save data
    temp.table.list <- lapply(split(temp.cgh.df, temp.cgh.df$Study.id), function(x){
                          t(table(factor(x$F.stat, levels=c("+", "-")), 
                                  factor(x$ci.stat, levels=c("High-CI", "Low-CI"))))
                        })
    cgh.contigency.list[[each.feat]] <- do.call("rbind", lapply(temp.table.list, as.vector))
    cgh.contigency.list[[each.feat]] <- cbind(cgh.contigency.list[[each.feat]],
                                              rep(each.feat, dim(cgh.contigency.list[[each.feat]])[1]),
                                              rownames(cgh.contigency.list[[each.feat]]))
  }
  
  if(each.feat %in% 'Met.stat'){
    temp.cgh.df <- cgh.df[which(!cgh.df[,'Met.stat'] %in% 'Unk'),] #Remove the Unknown met status samples
    
    # Generate Contigency table and save data
    temp.table.list <- lapply(split(temp.cgh.df, temp.cgh.df$Study.id), function(x){
                          t(table(factor(x$Met.stat, levels=c("Yes", "No")), 
                                  factor(x$ci.stat, levels=c("High-CI", "Low-CI"))))
                        })
    cgh.contigency.list[[each.feat]] <- do.call("rbind", lapply(temp.table.list, as.vector))
    cgh.contigency.list[[each.feat]] <- cbind(cgh.contigency.list[[each.feat]],
                                              rep(each.feat, dim(cgh.contigency.list[[each.feat]])[1]),
                                              rownames(cgh.contigency.list[[each.feat]]))
    
  }
}

#####  Format the cumulative contigency table dataframe
cgh.contigency.df <- do.call("rbind", cgh.contigency.list)
cgh.contigency.df <- data.frame(cgh.contigency.df, stringsAsFactors=FALSE, check.names=FALSE)
colnames(cgh.contigency.df) <- c("hipos", "lowpos", "hineg", "lowneg", "alloc", "author")
cgh.contigency.df$author <- gsub("\\.", " et al, ", cgh.contigency.df$author)
cgh.contigency.df$hipos <- as.integer(cgh.contigency.df$hipos)
cgh.contigency.df$hineg <- as.integer(cgh.contigency.df$hineg)
cgh.contigency.df$lowpos <- as.integer(cgh.contigency.df$lowpos)
cgh.contigency.df$lowneg <- as.integer(cgh.contigency.df$lowneg)


#####  Use the metafor package to plot the forestplot
pdf("cghForestplot.OR.post_review.pdf", width=10)
cex.val <- 0.65
meta.meas="OR"  # OR = log odds ratio, RR = log relative risk
#Generate the random model to accomodate for mutliple studies
res <- rma(ai=hipos, bi=hineg, ci=lowpos, di=lowneg, data=cgh.contigency.df, measure=meta.meas,
           slab=author, method="REML")
#Generate forest plots for the different features
par.mar <- par("mar")
op.par <- par()
split.screen(matrix(c(0, 0.8, 0, 1, 0.8, 1.0, 0, 1), byrow=TRUE, ncol=4))
screen(1)
par(mar=c(5.1, 4.1, 4.1, 0))
res$slab <- gsub("\\..*$", "", res$slab, perl=TRUE)
res$slab <- gsub("NETSEQ et al, ext", "NETSeq_Validation", res$slab, perl=TRUE)
res$slab <- gsub("NETSEQ", "NETSeq_Discovery", res$slab, perl=TRUE, ignore.case = FALSE)
# Generate row-indices for the 3 groups: F.stat, Met.stat, NET
table.idx <- table(cgh.contigency.df$alloc)
table.idx <- table.idx[unique(cgh.contigency.df$alloc)]
tbl.inf <- getTableIndices(table.idx)


forest(res, xlim=c(-16, 8), at=log(c(.05, .25, 1, 4, 20)), atransf=exp,
       ilab=cbind(cgh.contigency.df$hipos, cgh.contigency.df$hineg, cgh.contigency.df$lowpos, cgh.contigency.df$lowneg),
       ilab.xpos=c(-9.5,-8,-6,-4.5), cex=cex.val, ylim=c(-1, tbl.inf$ylim),
       rows=tbl.inf$rows,
       xlab=if(meta.meas %in% 'OR') "Odds Ratio" else "Relative Risk", mlab="RE Model for All Studies", psize=1)

### set font expansion factor (as in forest() above) and use bold italic
### font and save original settings in object 'op'
op <- par(cex=cex.val, font=4)

### add text for the subgroups
text(-16, tbl.inf$subs, pos=4, c("Metastasis Status",
                                 "Functional Status",
                                 "NET Type"))
### switch to bold font
par(font=2)

### add column headings to the plot
text(c(-9.5,-8,-6,-4.5), tbl.inf$cols['G3'], c("Met+", "Met-", "Met+", "Met-"))
text(c(-9.5,-8,-6,-4.5), tbl.inf$cols['G2'], c("F", "NF", "F", "NF"))
text(c(-9.5,-8,-6,-4.5), tbl.inf$cols['G1'], c("PNET", "GINET", "PNET", "GINET"))
text(c(-8.75,-5.25),     tbl.inf$cols['Title'], c("High-CI", "Low-CI"))
text(-16,                tbl.inf$cols['Title'], "Author(s) and Year",     pos=4)
text(8,                  tbl.inf$cols['Title'], paste(if(meta.meas %in% 'OR') "Odds Ratio" else "Relative Risk", "[95% CI]", sep=" "), pos=2)

### set par back to the original settings
par(op)

### fit random-effects model in the three subgroups
res.net <- rma(ai=hipos, bi=hineg, ci=lowpos, di=lowneg, data=cgh.contigency.df, measure=meta.meas,
               subset=(alloc=="NET"), method="REML")
res.fstat <- rma(ai=hipos, bi=hineg, ci=lowpos, di=lowneg, data=cgh.contigency.df, measure=meta.meas,
                 subset=(alloc=="F.stat"), method="REML")
res.met <- rma(ai=hipos, bi=hineg, ci=lowpos, di=lowneg, data=cgh.contigency.df, measure=meta.meas,
               subset=(alloc=="Met.stat"), method="REML")


### add summary polygons for the three subgroups
addpoly(res.net, row=tbl.inf$poly['G1'], cex=cex.val, atransf=exp, mlab="RE Model for Subgroup")
addpoly(res.fstat, row= tbl.inf$poly['G2'], cex=cex.val, atransf=exp, mlab="RE Model for Subgroup")
addpoly(res.met, row= tbl.inf$poly['G3'], cex=cex.val, atransf=exp, mlab="RE Model for Subgroup")


screen(2)
# Hardcoded extension to forestplot
x.mar <- par()$mar
par.mar.adj <- par.mar - c(0, 4.1, 1, 0)
par(mar = par.mar.adj)
plot(0, type='n', ylim=c(-1,tbl.inf$ylim), xlim=c(0,3), xlab="", ylab="", axes=FALSE)
abline(h = c(tbl.inf$ylim-2, 0), lty = 'solid', col = "black")

# Add in values for pval 
text(0.75, tbl.inf$poly['G1'], labels =specify_decimal(round(res.net$pval,3),3), cex=cex.val)
text(0.75, tbl.inf$poly['G2'], labels =specify_decimal(round(res.fstat$pval,3),3), cex=cex.val)
text(0.75, tbl.inf$poly['G3'], labels =specify_decimal(round(res.met$pval,3),3), cex=cex.val)
# Add in values for I^2
text(2, tbl.inf$poly['G1'], labels =specify_decimal(round(res.net$I2,2),2), cex=cex.val)
text(2, tbl.inf$poly['G2'], labels =specify_decimal(round(res.fstat$I2,2),2), cex=cex.val)
text(2, tbl.inf$poly['G3'], labels =specify_decimal(round(res.met$I2,2),2), cex=cex.val)

### set font expansion factor (as in forest() above) and use bold 
par(cex=cex.val, font=4)
par(font=2)
text(0.75, tbl.inf$poly['Title'], "pval")
text(2, tbl.inf$poly['Title'], "I^2")
close.screen(all=TRUE)
dev.off()

save(cgh.contigency.df, file="~/git/net-seq/cgh_analysis/cgh_plotter/output/cgh.contigency.post_review.Rdata")
 