# Uses equations described in NET-seq paper to calculate purity (pur.t) given an observed allelic fraction (Falt)
# and simulated ploidy in tumor (Pt), and theoretical fractions in tumor (Ft) and normal (Fn; where Fn=0.5 is germline)
calcPurity <- function(Ft=1, Pt=2, Fn=0.5, Pn=2, Falt=1){
  #pur.val <- Pn*(Fn - Falt) / ((Pt*Ft*(Falt - 1)) + (Pn*(Fn - Falt)))  # original incorrect formula
  #pur.t <- (Fn - Falt) / (Fn - Ft)  # Incorrect formula that doesnt consider ploidy
  #pur.t <- (2 * (Falt - Fn)) / ((Pt * Ft) - (2 * Fn))   # Incorrect same formula as above, but factors in ploidy
  
  Falt <- as.numeric(as.character(Falt))
  Pn <- as.integer(as.character(Pn))
  Pt <- as.integer(as.character(Pt))
  
  A = (Ft * (Pt / 2))
  a = (Fn * (Pn / 2))
  R = ((1-Ft) * (Pt / 2))
  r = ((1-Fn) * (Pn / 2))
  
  norm.val <- (a * (Falt - 1) + (Falt * r))
  total.val <- (a * (Falt - 1)) - (A * Falt) + (A) + (Falt * (r - R))
  pur.t <- norm.val / total.val
  
  return(pur.t)
}

# Uses equations described in NET-seq paper to calculate the expected allelic fraction (Falt) for a given 
# purity/ploidy and theoretical AF (Ft)
calcFalt <- function(Ft=1, Pt=2, Fn=0.5, Pn=2, pur.t=1){
  Pn <- as.integer(as.character(Pn))
  pur.t <- as.numeric(as.character(pur.t))
  Pt <- as.integer(as.character(Pt))
  
  A = (Ft * (Pt / 2)) * pur.t
  a = (Fn * (Pn / 2)) * (1 - pur.t)
  R = ((1-Ft) * (Pt / 2)) * pur.t
  r = ((1-Fn) * (Pn / 2)) * (1 - pur.t)
  
  Falt = (A + a) / (A + a + R + r)

  return(Falt)
}

# Uses equations described in NET-seq paper to calculate theoretical allelic fractions (Ft) for a given model and a
# corresponding observed allelic fraction (Falt)
calcTAF <- function(Falt=1, Pn=2, pur.t=1, Pt=2){
  Falt <- as.numeric(as.character(Falt))
  Pn <- as.integer(as.character(Pn))
  pur.t <- as.numeric(as.character(pur.t))
  Pt <- as.integer(as.character(Pt))
  
  Pt.pu <- Pt * pur.t
  
  Ft <- (Falt * ((-Pn*pur.t) + Pn + Pt.pu)) / Pt.pu
  
  return(Ft)
}

# Uses equations described in NET-seq paper to calculate Ploidy (each.cn) and B-alleles (Bt) given a known purity (pt)
# and theoretical values of normal allelic fraction (Fn; where Fn=0.5 is germline) and observed allelic fraction (Falt)
calcPloidy <- function(Falt=1, Fn=0.5, pt=1){
  Ft <- (Falt + Fn*(pt - 1))/pt
  cn.states <- c(1:10)
  viable.cn <- c()
  for(each.cn in cn.states){
    Bt <- (Ft * each.cn)
    if((Bt %% 1) == 0) viable.cn <- c(viable.cn, each.cn)
  }
  return(data.frame("copy.state"=viable.cn,
                    "B.val"=(viable.cn * Ft),
                    "Ft"=rep(Ft, length(viable.cn))))
}


# Simulates all possible CN-states and theoretical AF assuming 100% purity and runs each variant
# through the calcPurity() function
genDf <- function(af.list, id='null'){
  # Generates all possible AF for ploidy up to 10
  all.FN <- c("_somatic"=0, "_germline"=0.5)
  all.FN <- c("_somatic"=0)
  
  ploidy.models <- lapply(c(1:2), function(norm.ploidy){
    ploidy.models <- data.frame()
    for(i in c(1:8)){
      alt.allele <- seq(1,i)
      cn.stat <- rep(i, length(alt.allele))
      Ft <- round(alt.allele / cn.stat,3)
      ploidy.models <- rbind(ploidy.models, data.frame("Bt"=alt.allele,
                                                       "Pt"=cn.stat,
                                                       "Ft"=Ft))
    }
    ploidy.models$Pn <- rep(norm.ploidy, nrow(ploidy.models))
    return(ploidy.models)
  })
  ploidy.models <- do.call("rbind", ploidy.models)
  
  af.vals <- lapply(af.list, function(af.tum,FN){
    t(as.matrix(apply(ploidy.models, 1, function(each.model){
      round(calcPurity(Ft=each.model['Ft'], 
                       Pt=each.model['Pt'], 
                       Fn=FN, 
                       Pn=each.model['Pn'], 
                       Falt=af.tum), 3)
    })))
  }, FN=all.FN)
  
  # rename columns and rows and appends AF values
  af.vals <- lapply(names(af.vals), function(i) {
    af.vals[[i]] <- cbind(rep(round(af.list[[i]],3), length(all.FN)), af.vals[[i]])
    colnames(af.vals[[i]]) <- c("AF", paste(ploidy.models$Bt, ploidy.models$Pt, sep="-"))
    rownames(af.vals[[i]]) <- if(length(all.FN) > 1) paste0(i, names(all.FN)) else i
    af.vals[[i]]})
  names(af.vals) <- names(af.list)
  
  # Filter impossible purities
  af.vals <- lapply(af.vals, function(i){
    i[i < 0] <- NA
    i[i > 1] <- NA
    # na.idx <- apply(i, 2, function(j) all(is.na(j)))
    # if(length(which(na.idx)) > 0) i <- i[,-which(na.idx)]
    return(i)
  })
  
  #Tag on the sample ID and gene IDs
  af.vals <- as.data.frame(do.call("rbind", af.vals))
  af.vals <- cbind("ID"=rep(id, nrow(af.vals)),
                   "Gene"=rownames(af.vals),
                   af.vals)
  
  return(af.vals)
}