"""
## Date Created: 2-Dez-19
## Author: N. S. Araujo
## Figure of TEs types in the genome
"""

## Plot 1 - figure with the general proportions with data from RepeatMasker table output.

barplot(38.68, horiz = T, xlim=c(0,100),ylim=c(0,2), xaxt='n',
        las=1,cex.names=1, xlab="Percent of the genome", font.lab=2,col="white")
axis(side=1, at=seq(0,100, 10),font.axis=1)
barplot(1.31+2.32+7.17+11.23, horiz = T, add = T, axes = F,
        col="lightblue")
barplot(1.31+2.32+7.17, horiz = T, add = T, axes = F,
        col="lightsalmon")
barplot(1.31+2.32, horiz = T, add = T, axes = F,
        col="salmon3")
barplot(1.31, horiz = T, add = T, axes = F,
        col="grey")
legend("topright", fill = c("white","lightblue","lightsalmon","salmon3","grey"), legend = c("unclassified", "class II", "class I", "non-interspersed", "rolling-circles"), horiz = F, bty = "n") 


## Plot 2 - figure with all TEs family classification

## Al TEs in the genome
all_TEs <- read.csv("TDIV_genome_05_2019.fasta_formated.txt", header=T, sep = "\t",row.names=NULL, stringsAsFactors = FALSE)

## Filtering classification data
length(all_TEs$repeat_class.family)
TE_class.family <- as.character(all_TEs$repeat_class.family)
TE_class.family <- table(TE_class.family)
sum(TE_class.family)

# set colours according to TE category
col_TEs <- c("lightblue", "lightblue", "lightsalmon",
             "lightsalmon","lightblue","lightsalmon",
             "lightblue","lightblue","lightsalmon",
             "lightblue","lightsalmon", "black",
             "lightsalmon","lightblue","lightblue",
             "lightblue","lightblue","lightsalmon",
             "lightblue","lightblue","lightblue",
             "lightblue", "lightsalmon","lightsalmon",
             "lightblue","lightblue","lightsalmon",
             "lightsalmon","lightsalmon","lightsalmon",
             "lightblue","lightblue","lightsalmon",
             "lightsalmon","lightblue","lightblue",
             "lightsalmon","lightblue","lightblue",
             "lightblue","lightsalmon","lightsalmon",
             "lightblue","lightsalmon","salmon3",
             "salmon3","lightblue","lightsalmon",
             "lightblue","lightsalmon","lightsalmon",
             "lightblue","lightsalmon","lightsalmon",
             "lightsalmon","lightblue","lightsalmon",
             "salmon3","lightblue","grey",
             "lightblue","lightblue","salmon3",
             "lightblue","black")


## Save the plot in pdf

pdf.name <- "TEs_familiesINgenome.pdf"
pdf (pdf.name, width=12, height=5, compress = F)

par(mar=c(10,1.1,4.1,8.1))
plot(sort(TE_class.family), las=2, xlab="", ylab="", 
     yaxt='n',xaxt='n', ylim = c(0,360000))
axis(side=4, at=seq(0,360000, 50000),
     font.axis=1, las=2, cex.axis=0.7)
axis(side=1, at=seq(1,65, 1), 
     labels=c(rep('',65)))
for (j in 1:65){
  axis(side=1, at=j, col.axis=col_TEs[j], names(sort(TE_class.family))[j], las=2,
       font.axis=2,cex.axis=0.6) # Add habitat as labels, each with corresponding color, on the left margin
}

mtext("# of annotated entries", side=4, line=3, cex.lab=1,las=3)
mtext("repeat class/family", side=1, line=7, cex.lab=1,las=1)


dev.off()

