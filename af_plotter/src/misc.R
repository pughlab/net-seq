fixFactors <- function(x.df, col.n, int.or.num='integer'){
  if(int.or.num %in% 'integer'){
    for(e.col in col.n){
      x.df[,e.col] <- as.integer(as.character(x.df[,e.col]))
    }
  } else if (int.or.num %in% 'numeric'){
    for(e.col in col.n){
      x.df[,e.col] <- as.numeric(as.character(x.df[,e.col]))
    }
  }
  return(x.df)
}

raw.t.test <- function(ds1, ds2){
  if(!is.na(ds1['x']) & !is.na(ds2['x'])){
    if(ds1['x'] > ds2['x']){
      ds3 <- ds1
      ds1 <- ds2
      ds2 <- ds3
    }
    t.val <- (ds1['x']- ds2['x']) / sqrt(ds1['sd']^2/ds1['n']+ds2['sd']^2/ds2['n'])
    A=ds1['sd']^2/ds1['n']
    B=ds2['sd']^2/ds2['n']
    df=(A+B)^2/(A^2/(ds1['n']-1)+B^2/(ds2['n']-1))
    pval <- pt(t.val, df) * 2
  } else {
    pval <- NA
  }
  return(pval)
}
