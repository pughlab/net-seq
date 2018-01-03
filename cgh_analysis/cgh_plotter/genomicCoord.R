#### Get's the base-pair positions for a string of karyotype locations
getCytoPos <- function(cn.string, hg19.cytoband){
  cn.string <- gsub("[\" ]", "", cn.string)
  cn.string <- gsub("None", "", cn.string)
  cn.v <- unlist(strsplit(cn.string, split = ","))
  cn.list <- lapply(cn.v, function(x) splitCytoLoc(x))
  cn.coord.list <- lapply(cn.list, function(x){chr.sub <- hg19.cytoband[grep(paste(x['chr'], "$", sep=""), hg19.cytoband$chrom),]
                                               spos <- formatCytoPos(x['start'], 'start', chr.sub)
                                               epos <- formatCytoPos(x['end'], 'end', chr.sub)
                                               i <- c("chr" = x["chr"],
                                                      "start" = spos,
                                                      "end" = epos)
                                               return(i)
                                               })
  return(cn.coord.list)
}

#### Gets the genomic coordinates for a cytoband locations (e.g. chr1q21.3)
formatCytoPos <- function(pos, ord, sub.df){
 spos <- NA   # Will store the genomic location for a given cytoband
 
 #If the entire chromosome is altered and no start/end is given
 if(nchar(pos) == 0){
   if(ord == "start"){
     spos <- sub.df[1, 'chromStart']
   } else if (ord=='end') {
     spos <- sub.df[dim(sub.df)[1], 'chromEnd']
   } 
   
 # If the entire arm is altered or an actual position is given
 } else if(length(grep("[pq]", pos))>0) {
   # Tests if the qterminal "qter" 
   if(length(grep("qter", pos)) > 0){
     pos <- gsub("ter", "", pos)
     ord <- "end"
   }
   # or the pterminal "pter" is given
   if(length(grep("pter", pos)) > 0){
     pos <- gsub("ter", "", pos)
     ord <- 'start'
   }
   
   # Checks if the position is found in the cytoband locations
   fix.cnt <- 1
   while(length(grep(paste("^", pos, sep=""), sub.df$name)) == 0){
     print(paste("Position not found: ", sub.df[1, 'chrom'], ": ", pos, sep=""))
     if(fix.cnt == 1){
       new.pos <- gsub("\\.[0-9]+$", "" , pos)
     } else if(fix.cnt == 2){
       new.pos <- gsub("[0-9]+$", "", pos)
     } else {
       warning(paste("Could not find a proper position for ", pos, sep=""))
     }
     
     print(paste("Converting ", pos, " to ", new.pos, sep=""))
     pos <- new.pos
     fix.cnt <- fix.cnt + 1
   }
   
   sub.df.rows <- grep(paste("^", pos, sep=""), sub.df$name)
   if(ord == "start"){
     spos <- sub.df[min(sub.df.rows), 'chromStart']
   } else if (ord=='end') {
     spos <- sub.df[max(sub.df.rows), 'chromEnd']
   } else {
     warning("No start/end order given.")
   }
   
 #If the centrosome is the starting/ending point
 } else if(pos == 'cen'){
   pos <- "acen"
   sub.df.rows <- which(sub.df$gieStain %in% pos)
   if(ord == "start"){
     spos <- sub.df[min(sub.df.rows), 'chromEnd']
   } else if (ord=='end') {
     spos <- sub.df[max(sub.df.rows), 'chromStart']
   } else {
     warning("No start/end order given.")
   }
 
 #If no criteria was matched  
 } else {
   warning("Pos did not match any of the given options")
 }
 
 return(spos)
}


#### Splits a string from 20q12-q13.1 to chr20; q12; q13.1
splitCytoLoc <- function(loc){
  # -------debug code
  #   iter <- c("20q12-q13.1",
  #             "5cen-q23",
  #             "7p",
  #             "8",
  #             "20q13.3")
  #   for(loc in iter){
  #     
  #   }
  #Regex for the Chr num
  chr <- gsub("[pq(cen)].*?$", "", loc)
  chr <- gsub("^chr", "", chr)
  chr.f <- paste("chr", chr, sep="")
  
  #Regex for the start position
  start <- gsub("(-.+)?", "", loc)
  start <- gsub(paste("^", chr, sep=""), "", start)
  
  #Regex for the start position
  if(length(grep("-", loc) ) > 0){  
    end <- gsub(".+?-", "", loc)
  } else {  # No end position given, e.g. "8p"
    end <- start
  }
  
  x <- c("chr"=chr.f, "start"=start, "end"=end)
  return(x)
}

