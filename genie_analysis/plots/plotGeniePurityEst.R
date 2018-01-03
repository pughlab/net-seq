purity.df <- read.table("~/git/net-seq/genie_analysis/data/theoretical_estimated_purity_frac.txt",
                        sep="\t", stringsAsFactors = FALSE, check.names = FALSE, header = TRUE)
dir.create("~/Desktop/netseq/genie", showWarnings = FALSE, recursive = TRUE)
pdf("~/Desktop/netseq/genie/purity_est.noLabel.pdf", width=5, height=5)
purity.df[is.na(purity.df)] <- 0
with(purity.df, plot(Purity_frac, denovo_purity_est,
                     xlim=c(0,1.15), ylim=c(0,1.15),
                     ylab="Estimated Purity",
                     xlab="Pathologist Purity",
                     main="Denovo purity estimates"))
lines(x=c(-1,2), y=c(-1,2), lty=3)
 # with(purity.df, text(x=(Purity_frac + 0.05), 
 #                      y=denovo_purity_est, 
 #                      labels = gsub("GENIE_Pnet", "Pnet", UID), adj=0, cex=0.5))
 
purity.df[is.na(purity.df)] <- 0
with(purity.df, plot(Purity_frac, manul_purity_est,
                    xlim=c(0,1.15), ylim=c(0,1.15),
                    ylab="Estimated Purity",
                    xlab="Pathologist Purity",
                    main="Pathologist-guided purity estimates"))
lines(x=c(-1,2), y=c(-1,2), lty=3)
# with(purity.df, text(x=(Purity_frac + 0.05), 
#                      y=manul_purity_est, 
#                      labels = gsub("GENIE_Pnet", "Pnet", UID), adj=0, cex=0.5))
dev.off()
