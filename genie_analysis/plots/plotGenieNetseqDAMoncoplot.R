library(ggplot2)
library(reshape)
library(scales)
col.idx <- list("Ins"="#82CABC",
                "Del"="#F26B67",
                "Non"="#B1ADCE",
                "Mis"="#72A2C6",
                "Dtrunc"="#000000",
                "Ptrunc"="#7B8390",
                "Inframe"="#FF9C2D",
                "0"="gray92")

legend.df <- data.frame("SampleID"=paste0("legend", c(1:length(col.idx))),
                        "MEN1"=names(col.idx),
                        "DAXX"=rep(0, length(col.idx)),
                        "ATRX"=rep(0, length(col.idx)))

net.df <- data.frame("SampleID"=c("NET-001-T_2", "NET-003-T_1", "NET-003-T_2", "NET-008-T_2", "NET-009-T_1", "NET-009-T_2"),
                     "MEN1"=c("Del", "Del", "Del", "Ins", "Non", "Non"),
                     "DAXX"=c("Mis", "Mis", 0, 0, 0, 0),
                     "ATRX"=c(0, 0, 0, "Ins", "Mis", "Mis"))
dam.df <- data.frame("SampleID"=paste0("GENIE_Pnet_", c(9:37)),
                     "MEN1"=c(rep("Dtrunc", 2), "Mis", rep("Dtrunc", 2), "Mis", rep("Dtrunc", 2), 
                              "Inframe", "Dtrunc", 0, 0, rep("Dtrunc", 3), "Mis", rep("Dtrunc", 5), 
                              0, rep("Dtrunc", 3), 0, rep("Dtrunc", 3)),
                     "DAXX"=c(0,0,"Mis", "Ptrunc", 0, rep("Ptrunc", 3), 0, 0, "Ptrunc", 0, 0, 0, 0, 
                              "Ptrunc", 0, "Mis", rep("Ptrunc", 2), 0, 0, rep("Ptrunc", 4), 0, "Mis", 0),
                     "ATRX"=c(rep(0,4), "Mis", rep(0,3), "Inframe", 0,0, rep("Mis", 2), "Dtrunc", "Mis",
                              0,0,"Mis", 0,0, rep("Dtrunc", 2), rep(0,4), "Dtrunc", 0,0))
nondam.df <- data.frame("SampleID"=paste0("GENIE_Pnet_", c(40:55)),
                     "MEN1"=rep(0,16),
                     "DAXX"=rep(0,16),
                     "ATRX"=rep(0,16))


plotOncoMap(legend.df, h=8, w=4, file.id="~/Desktop/legend.pdf")
plotOncoMap(net.df, h=4, w=4, file.id="~/Desktop/net.pdf")
plotOncoMap(dam.df, h=19, w=4, file.id="~/Desktop/dam.pdf")
plotOncoMap(nondam.df, h=16, w=4, file.id="~/Desktop/nondam.pdf")

plotOncoMap <- function(x.df, h, w, file.id){
  pdf(file.id, height=h, width=w)
  split.screen(c(1,3))
  screen(1) #MEN1
  plotRect(x.df, 'MEN1')
  
  screen(2) #DAXX
  plotRect(x.df, 'DAXX')
  
  screen(3) #ATRX
  plotRect(x.df, 'ATRX')
  close.screen(all.screens=TRUE)
  dev.off()
}




plotRect <- function(x.df, col.id){
  row.idx <- split.screen(c(nrow(x.df), 1))
  for(each.row in c(1:nrow(x.df))){
    screen(row.idx[each.row])
    par(mar=c(0.2, 1, 0.2, 1))
    plot(0, type='n', xlim=c(0,1), ylim=c(0,1), axes = FALSE, ylab='', xlab='', xaxt='n', yaxt='n')
    rect(0, 0, 1, 1, col = col.idx[[as.character(x.df[each.row, col.id])]], border = NA)
  }
}

pdf("~/Desktop/gain_loss.pdf")
colfunc <- colorRampPalette(c("blue", "white", "red"))
legend_image <- as.raster(matrix(colfunc(100), ncol=1))
plot(x = c(0,2), y = c(0,1),
     type = 'n', axes = F,xlab = '', ylab = '', main = 'legend title')
rasterImage(legend_image, 0, 0, 1, 1)
dev.off()