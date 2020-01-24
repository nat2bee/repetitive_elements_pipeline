"""
## Date Created: 2-Dez-19
## Author: N. S. Araujo
## Position of TEs around genes including all genes and only the DET in the same plot
"""


### Run 'TEs_around_gene.R' to get TEs around genes in the two sets of genes:

## 1- gene DE in diapausa
DET_data2plot <- data2plot
DET_table_data2plot <- table(data2plot)
DET_TEinGENE <- TEinGENE

## 2- all genes in the genome
allgenes_data2plot <- data2plot
allgenes_table_data2plot <- table(data2plot)
allgenes_TEinGENE <- TEinGENE


## Plot the results
pdf.name <- "TEs_close2_genes.pdf"
pdf (pdf.name, width=14, height=6, compress = F)
par(mfrow=c(1,2),mar=c(5.1,1,5.1,1))

## all genomic genes
hist(allgenes_data2plot,  breaks=1000, probability = T, xlim = c(-60, 60), ylim = c(0.05,0.15),
     border=rgb(0,0,0,0.5), main="", xlab="",xaxt='n', yaxt='n')
axis(side=1, at=seq(-60,60, 10), labels = c(seq(-6000,6000,1000)),
     cex.axis=0.7) 
mtext("position upstream and downstream of the genes (bp)", side=1, line=3, cex.lab=1,las=1)
mtext("all genes", side=3, line=1, cex.lab=1,las=1)

## DET in diapause
hist(DET_data2plot,  breaks=1000, probability = T, xlim = c(-60, 60), ylim = c(0.05,0.15),
     border=rgb(1,0,0,0.5), main= "", xlab="",xaxt='n', yaxt='n')
axis(side=1, at=seq(-60,60, 10), labels = c(seq(-6000,6000,1000)),
     cex.axis=0.7) 
mtext("position upstream and downstream of the genes (bp)", side=1, line=3, cex.lab=1,las=1)
mtext("differentially expressed in diapause", side=3, line=1, cex.lab=1,las=1)

## central axis
axis(side=2, at=seq(0.05,0.15, 0.025),
     font.axis=1, las=2, cex.axis=0.7)

dev.off()
