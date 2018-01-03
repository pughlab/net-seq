#### Draws Losses and Gains rectangles for each sample on an existing plot
plotCNA <- function(cna.info, cna, cna.col="grey", ucsc.chrom, t.genome.size){
  #cna.info contains a list of genomic coordinates as a labelled vector in each element
  plotRect <- function(x, ucsc.chrom, t.genome.size){
    if(!x['chr.chr'] == 'chrX' & !x['chr.chr'] == 'chrY'){
      chr.factor <- ucsc.chrom[which(ucsc.chrom$chrom == x['chr.chr']), 'chr.st.size']
      
      start.coord <- -as.integer(x['start']) - chr.factor
      end.coord <- -as.integer(x['end']) - chr.factor
      
      rect(xleft = 0, ybottom = end.coord, xright = 1, ytop = start.coord, col = cna.col)
    }
  }
  
  if(length(cna.info[[1]]) > 0){ 
    lapply(cna.info[[1]], function(x) plotRect(x, ucsc.chrom, t.genome.size))
  }
}

#### Converts a given string into a set of colours based on the levels/factors
convertFactorToColor <- function(group.v, cgh.pal, ret.ord=0){
  col.palettes <- getCghColPalette()
  group.v <- factor(group.v)
  col.ord <- match(group.v, levels(group.v))
  
  col.palette <- col.palettes[[cgh.pal]]
  #col.palette <- pal(length(levels(group.v)))
  col.ord <- col.palette[col.ord]
  if(ret.ord == 0){
    return(col.ord)
  } else if (ret.ord == 1){
    uniq.df <- unique(data.frame(col.ord, group.v))
    col.list <- list("color"=as.character(uniq.df$col.ord),
                     "id"=as.character(uniq.df$group.v))
    return(col.list)
  }
  
}


#### Colour Palettes for different groups
getCghColPalette <- function(){
  require(gplots)
  col.palettes <<- list()
  # # NET
  # > levels(factor(cgh.df$NET))
  # [1] "GINET" "LNET"  "PNET" 
  col.palettes[['NET']] <- alpha(c(gray(0.9), redgreen(1), gray(0)), 0.75)
  
  # # Study.id
  # > levels(factor(cgh.df$Study.id))
  # [1] "Floridia.2005" "Haugvik.2014"  "Nagano.2007"   "NETSEQ"        "NETSEQ.ext"    "Speel.2001"    
  # [7] "Stumpf.2000"   "Terris.1998"  "Tonnies.2000"  "Zhao.2001" 
  col.palettes[['Study.id']] <- alpha(rainbow(10), 0.75)
  
  # # F.stat
  # > levels(factor(cgh.df$F.stat))
  # [1] "-"   "+"   "Unk"
  col.palettes[['F.stat']] <- alpha(c(gray(0.9), gray(0), redgreen(1)), 0.75)
  
  # # Met.stat
  # > levels(factor(cgh.df$Met.stat))
  # [1] "No"  "Unk" "Yes"
  col.palettes[['Met.stat']] <- alpha(c(gray(0.9), redgreen(1), gray(0)), 0.75)
  
  # # ci.stat
  # > levels(factor(cgh.df$ci.stat))
  # [1] "High-CI" "Low-CI" 
  col.palettes[['ci.stat']] <- alpha(gray(c(0,0.9)), 0.75)
  
  return(col.palettes)
}
