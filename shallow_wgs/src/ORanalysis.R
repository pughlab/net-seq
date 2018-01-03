specify_decimal <- function(x, k) format(round(x, k), nsmall=k)

genORTable <- function(each.group){
  group.id <- each.group[['group.id']]
  group <- each.group[['group']]
  
  colnames(snp.deficient.df) <- gsub("^X", "NET-2-", colnames(snp.deficient.df))
  colnames(cn.state.df) <- gsub("^X", "NET-2-", colnames(cn.state.df))
  
  group.list <- split(anno.df, anno.df[,group])
  group.stats <- lapply(group.id, function(each.group.id){
    # Isolate each groups samples
    print(paste("Working with: ", each.group.id, sep=""))
    group.col <- which(colnames(snp.deficient.df) %in% group.list[[each.group.id]][,'sample_id'])
    group.deficient.df <- snp.deficient.df[,group.col]
    group.cn.df <- cn.state.df[,group.col]
    
    # Count how many Het-SNP deficient/enriched chromosomes there are in each group
    def.count <- apply(group.deficient.df, 1, function(x) sum(x, na.rm=TRUE))
    nonDef.count <- apply(group.deficient.df, 1, function(x) sum(abs(x - 1), na.rm=TRUE))
    
    #Obtain CN state tally for each group:    
    if(dim(group.cn.df)[2] == 0){
      def.mean.sd <- rep(0, nrow(group.cn.df))
      nonDef.mean.sd <- rep(0, nrow(group.cn.df))
    } else {
      def.mean.sd <- apply(group.cn.df[,which(group.deficient.df[each.chr,] == 1), drop=FALSE], 
                           1, getMeanAndSd)
      nonDef.mean.sd <- apply(group.cn.df[,which(group.deficient.df[each.chr,] == 0), drop=FALSE],
                              1, getMeanAndSd)
    }
    
    
    
    group.deficient.df <- cbind(group.deficient.df, 
                                data.frame(snp.deficient=def.count,
                                           snp.enriched=nonDef.count),
                                data.frame(cn.deficient=def.mean.sd,
                                           cn.enriched=nonDef.mean.sd))
    return(list(af=group.deficient.df[,c('snp.deficient', 'snp.enriched')],
                cn=group.deficient.df[,c('cn.deficient', 'cn.enriched')]))
  })
  
  #   group.stats.def <- as.data.frame(cbind(group.stats[[1]]$af, group.stats[[2]]$af,
  #                                          group.stats[[1]]$cn, group.stats[[2]]$cn))
  #   colnames(group.stats.def) <- c("a.deficient", "a.enriched", 
  #                                  "b.deficient", "b.enriched",
  #                                  "a.def.cn", "a.enr.cn",
  #                                  "b.def.cn", "b.enr.cn")
  
  group.stats.def <- as.data.frame(cbind(group.stats[[1]]$af, group.stats[[2]]$af))
  colnames(group.stats.def) <- c("a.deficient", "a.enriched", 
                                 "b.deficient", "b.enriched")
  
  group.stats.def <- cbind(group.stats.def, 
                           alloc=rep(group), 
                           chrom=paste("chr", rownames(group.stats.def), sep=""))
  group.stats.def <- group.stats.def[c(dim(group.stats.def)[1]:1),]
  res <- rma(ai=a.deficient, bi=a.enriched, ci=b.deficient, di=b.enriched, 
             data=group.stats.def, measure="OR", slab=chrom, method="REML")
  
  
  pdf(paste("/Users/rquevedo/Desktop/PughLab/NET-seq/otb_samples/net_OTB/plots/chrOR.", paste(group, collapse="."), ".pdf", sep=""))
  
  split.screen(matrix(c(0, 0.8, 0, 1, 0.8, 1.0, 0, 1), byrow=TRUE, ncol=4))
  screen(1)
  par(mar=c(5.1, 4.1, 4.1, 0))
  forest(res, xlim=c(-12, 8), at=log(c(.05, .25, 1, 4, 20)), atransf=exp,
         ilab=cbind(group.stats.def$a.deficient, group.stats.def$a.enriched, 
                    group.stats.def$b.deficient, group.stats.def$b.enriched),
         #ilab.xpos=c(-22.5, -21, -18, -16.5, -13.5,-12,-9,-7.5), cex=.75, ylim=c(-1, 29),
         ilab.xpos=c(-9.5,-8,-6.5,-5), cex=.65, ylim=c(-1, 29),
         rows=c(3:24),
         xlab="Odds Ratio", mlab="Genome", psize=1)
  
  
  
  ### add text for the subgroups
  op <- par(cex=.65, font=4)
  text(-12, c(25), pos=4, c("Chrom"))
  
  ### add column headings to the plot
  par(font=2)
  text(c(-9.5,-8,-6.5,-5), 26, c("LOH", "Het", "LOH", "Het"))
  text(c(-8.75,-5.75),     28, group.id)
  text(8,                  28, "Odds Ratio [95% CI]", pos=2)
  
  screen(2)
  # Hardcoded extension to forestplot
  par(mar=c(5.1, 0, 3.1, 2.1))
  plot(0, type='n', ylim=c(-1,29), xlim=c(0,20), xlab="", ylab="", axes=FALSE)
  abline(h = c(29 - 2, 0), lty = 'solid', col = "black")
  text(10, -1, labels =specify_decimal(round(res$pval,3),3), cex=0.65)  
  op <- par(cex=.65, font=2)
  text(10,     28, 'p.value')
  
  close.screen(all.screens=TRUE)
  ### set par back to the original settings
  dev.off()
  par(op)
}