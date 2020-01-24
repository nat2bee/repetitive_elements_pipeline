"""
## Date Created: 26-Nov-19
## Author: N. S. Araujo
## Script used to filter repetitive elements around genes in the genome
"""

## Get a list of genes of interest and their position in the genome
DET_diapausa <- read.csv("Gene_position_DEG_diapause_Tdiv.txt", header=T, sep = "\t")
  
## Default TEs list output "TDIV_genome_05_2019.fasta.out" needs to be formatted to be used.
## In a text editor I've removed the header and then with a bash script I've formatted it.
"""
$ tr -s '[:blank:]' '|' < TDIV_genome_05_2019.fasta_header.out | sed -e 's/^|//g' | sed -e 's/|$//g' | sed -e 's/*$//g' | sed -e 's/|$//g' | tr -s '|' '\t' > TDIV_genome_05_2019.fasta_formated.txt
"""
## This is the input file with all repetitive elements position in the genome
## Openning it
all_TEs <- read.csv("TDIV_genome_05_2019.fasta_formated.txt", header=T, sep = "\t",row.names=NULL, stringsAsFactors = FALSE)

## getting position of the TEs in relation to gene (5kb window)
data2plot <- vector()
checked_references <- vector()
TEinGENE <- data.frame(matrix(ncol = as.numeric(length(colnames(all_TEs)))+1, nrow = 0))
colnames(TEinGENE) <- c(colnames(all_TEs), "Diapause_gene")
a <- 0


for(i in 1:length(DET_diapausa$Scaffold)){
  reference <- as.character(DET_diapausa$Scaffold[i])
  # Check if this scafold wasn't already checked 
  if(!(reference %in% checked_references)){
    checked_references <- c(checked_references,reference)
    gene_start <- DET_diapausa$Gene_start[i]
    gene_end <- DET_diapausa$Gene_end[i]
    # check the gene sense and invert if the oposite frame (just for graphic use)
    if(gene_start > gene_end){
      gene_end <- DET_diapausa$Gene_start[i]
      gene_start <- DET_diapausa$Gene_end[i]
    }
    TEsInRef <- all_TEs[all_TEs$query_sequence == reference,] # Filter TEs in the same scaffold
    # for all TEs in the same scaffold check if they occur in the gene region
    for(j in 1:length(TEsInRef$SW_score)){
      plot_pos <- FALSE
      TE_start <- TEsInRef$query_begin[j]
      TE_end <- TEsInRef$query_end[j]
      # excatly within the gene
      if((gene_start <= TE_start && TE_start <= gene_end) || (gene_start <= TE_end && TE_end <= gene_end)){
        plot_pos <- 0
      }
      
      # upstream the gene, up to 5kb
      else if((gene_end <= TE_start && (gene_end+5000) >= TE_start) || (TE_end >= gene_end && TE_end <= (gene_end+5000))){
        if((gene_end <= TE_start && (gene_end+100) >= TE_start) || (TE_end >= gene_end && TE_end <= (gene_end+100))) {
          plot_pos <- 1
        }
        else if((gene_end <= TE_start && (gene_end+200) >= TE_start) || (TE_end >= gene_end && TE_end <= (gene_end+200))){
          plot_pos <- 2
        }
        else if((gene_end <= TE_start && (gene_end+300) >= TE_start) || (TE_end >= gene_end && TE_end <= (gene_end+300))) {
          plot_pos <- 3
        }
        else if((gene_end <= TE_start && (gene_end+400) >= TE_start) || (TE_end >= gene_end && TE_end <= (gene_end+400))){
          plot_pos <- 4
        }
        else if((gene_end <= TE_start && (gene_end+500) >= TE_start) || (TE_end >= gene_end && TE_end <= (gene_end+500))) {
          plot_pos <- 5
        }
        else if((gene_end <= TE_start && (gene_end+600) >= TE_start) || (TE_end >= gene_end && TE_end <= (gene_end+600))){
          plot_pos <- 6
        }
        else if((gene_end <= TE_start && (gene_end+700) >= TE_start) || (TE_end >= gene_end && TE_end <= (gene_end+700))) {
          plot_pos <- 7
        }
        else if((gene_end <= TE_start && (gene_end+800) >= TE_start) || (TE_end >= gene_end && TE_end <= (gene_end+800))){
          plot_pos <- 8
        }
        else if((gene_end <= TE_start && (gene_end+900) >= TE_start) || (TE_end >= gene_end && TE_end <= (gene_end+900))) {
          plot_pos <- 9
        }
        else if((gene_end <= TE_start && (gene_end+1000) >= TE_start) || (TE_end >= gene_end && TE_end <= (gene_end+1000))){
          plot_pos <- 10
        }
        else if((gene_end <= TE_start && (gene_end+1100) >= TE_start) || (TE_end >= gene_end && TE_end <= (gene_end+1100))) {
          plot_pos <- 11
        }
        else if((gene_end <= TE_start && (gene_end+1200) >= TE_start) || (TE_end >= gene_end && TE_end <= (gene_end+1200))){
          plot_pos <- 12
        }
        else if((gene_end <= TE_start && (gene_end+1300) >= TE_start) || (TE_end >= gene_end && TE_end <= (gene_end+1300))) {
          plot_pos <- 13
        }
        else if((gene_end <= TE_start && (gene_end+1400) >= TE_start) || (TE_end >= gene_end && TE_end <= (gene_end+1400))){
          plot_pos <- 14
        }
        else if((gene_end <= TE_start && (gene_end+1500) >= TE_start) || (TE_end >= gene_end && TE_end <= (gene_end+1500))) {
          plot_pos <- 15
        }
        else if((gene_end <= TE_start && (gene_end+1600) >= TE_start) || (TE_end >= gene_end && TE_end <= (gene_end+1600))){
          plot_pos <- 16
        }
        else if((gene_end <= TE_start && (gene_end+1700) >= TE_start) || (TE_end >= gene_end && TE_end <= (gene_end+1700))) {
          plot_pos <- 17
        }
        else if((gene_end <= TE_start && (gene_end+1800) >= TE_start) || (TE_end >= gene_end && TE_end <= (gene_end+1800))){
          plot_pos <- 18
        }
        else if((gene_end <= TE_start && (gene_end+1900) >= TE_start) || (TE_end >= gene_end && TE_end <= (gene_end+1900))) {
          plot_pos <- 19
        }
        else if((gene_end <= TE_start && (gene_end+2000) >= TE_start) || (TE_end >= gene_end && TE_end <= (gene_end+2000))){
          plot_pos <- 20
        }
        else if((gene_end <= TE_start && (gene_end+2100) >= TE_start) || (TE_end >= gene_end && TE_end <= (gene_end+2100))) {
          plot_pos <- 21
        }
        else if((gene_end <= TE_start && (gene_end+2200) >= TE_start) || (TE_end >= gene_end && TE_end <= (gene_end+2200))){
          plot_pos <- 22
        }
        else if((gene_end <= TE_start && (gene_end+2300) >= TE_start) || (TE_end >= gene_end && TE_end <= (gene_end+2300))) {
          plot_pos <- 23
        }
        else if((gene_end <= TE_start && (gene_end+2400) >= TE_start) || (TE_end >= gene_end && TE_end <= (gene_end+2400))){
          plot_pos <- 24
        }
        else if((gene_end <= TE_start && (gene_end+2500) >= TE_start) || (TE_end >= gene_end && TE_end <= (gene_end+2500))) {
          plot_pos <- 25
        }
        else if((gene_end <= TE_start && (gene_end+2600) >= TE_start) || (TE_end >= gene_end && TE_end <= (gene_end+2600))){
          plot_pos <- 26
        }
        else if((gene_end <= TE_start && (gene_end+2700) >= TE_start) || (TE_end >= gene_end && TE_end <= (gene_end+2700))) {
          plot_pos <- 27
        }
        else if((gene_end <= TE_start && (gene_end+2800) >= TE_start) || (TE_end >= gene_end && TE_end <= (gene_end+2800))){
          plot_pos <- 28
        }
        else if((gene_end <= TE_start && (gene_end+2900) >= TE_start) || (TE_end >= gene_end && TE_end <= (gene_end+2900))) {
          plot_pos <- 29
        }
        else if((gene_end <= TE_start && (gene_end+3000) >= TE_start) || (TE_end >= gene_end && TE_end <= (gene_end+3000))){
          plot_pos <- 30
        }
        else if((gene_end <= TE_start && (gene_end+3100) >= TE_start) || (TE_end >= gene_end && TE_end <= (gene_end+3100))) {
          plot_pos <- 31
        }
        else if((gene_end <= TE_start && (gene_end+3200) >= TE_start) || (TE_end >= gene_end && TE_end <= (gene_end+3200))){
          plot_pos <- 32
        }
        else if((gene_end <= TE_start && (gene_end+3300) >= TE_start) || (TE_end >= gene_end && TE_end <= (gene_end+3300))) {
          plot_pos <- 33
        }
        else if((gene_end <= TE_start && (gene_end+3400) >= TE_start) || (TE_end >= gene_end && TE_end <= (gene_end+3400))){
          plot_pos <- 34
        }
        else if((gene_end <= TE_start && (gene_end+3500) >= TE_start) || (TE_end >= gene_end && TE_end <= (gene_end+3500))) {
          plot_pos <- 35
        }
        else if((gene_end <= TE_start && (gene_end+3600) >= TE_start) || (TE_end >= gene_end && TE_end <= (gene_end+3600))){
          plot_pos <- 36
        }
        else if((gene_end <= TE_start && (gene_end+3700) >= TE_start) || (TE_end >= gene_end && TE_end <= (gene_end+3700))) {
          plot_pos <- 37
        }
        else if((gene_end <= TE_start && (gene_end+3800) >= TE_start) || (TE_end >= gene_end && TE_end <= (gene_end+3800))){
          plot_pos <- 38
        }
        else if((gene_end <= TE_start && (gene_end+3900) >= TE_start) || (TE_end >= gene_end && TE_end <= (gene_end+3900))) {
          plot_pos <- 39
        }
        else if((gene_end <= TE_start && (gene_end+4000) >= TE_start) || (TE_end >= gene_end && TE_end <= (gene_end+4000))){
          plot_pos <- 40
        }
        else if((gene_end <= TE_start && (gene_end+4100) >= TE_start) || (TE_end >= gene_end && TE_end <= (gene_end+4100))) {
          plot_pos <- 41
        }
        else if((gene_end <= TE_start && (gene_end+4200) >= TE_start) || (TE_end >= gene_end && TE_end <= (gene_end+4200))){
          plot_pos <- 42
        }
        else if((gene_end <= TE_start && (gene_end+4300) >= TE_start) || (TE_end >= gene_end && TE_end <= (gene_end+4300))) {
          plot_pos <- 43
        }
        else if((gene_end <= TE_start && (gene_end+4400) >= TE_start) || (TE_end >= gene_end && TE_end <= (gene_end+4400))){
          plot_pos <- 44
        }
        else if((gene_end <= TE_start && (gene_end+4500) >= TE_start) || (TE_end >= gene_end && TE_end <= (gene_end+4500))) {
          plot_pos <- 45
        }
        else if((gene_end <= TE_start && (gene_end+4600) >= TE_start) || (TE_end >= gene_end && TE_end <= (gene_end+4600))){
          plot_pos <- 46
        }
        else if((gene_end <= TE_start && (gene_end+4700) >= TE_start) || (TE_end >= gene_end && TE_end <= (gene_end+4700))) {
          plot_pos <- 47
        }
        else if((gene_end <= TE_start && (gene_end+4800) >= TE_start) || (TE_end >= gene_end && TE_end <= (gene_end+4800))){
          plot_pos <- 48
        }
        else if((gene_end <= TE_start && (gene_end+4900) >= TE_start) || (TE_end >= gene_end && TE_end <= (gene_end+4900))) {
          plot_pos <- 49
        }
        else if((gene_end <= TE_start && (gene_end+5000) >= TE_start) || (TE_end >= gene_end && TE_end <= (gene_end+5000))){
          plot_pos <- 50
        }
        
      }
      
      # downstream the gene, up to 5kb
      else if((gene_start >= TE_start && (gene_start-5000) <= TE_start) || (TE_end <= gene_start && TE_end >= (gene_start-5000))){
        if((gene_start >= TE_start && (gene_start-100) <= TE_start) || (TE_end <= gene_start && TE_end >= (gene_start-100))) {
          plot_pos <- -1
        }
        else if((gene_start >= TE_start && (gene_start-200) <= TE_start) || (TE_end <= gene_start && TE_end >= (gene_start-200))){
          plot_pos <- -2
        }
        else if((gene_start >= TE_start && (gene_start-300) <= TE_start) || (TE_end <= gene_start && TE_end >= (gene_start-300))) {
          plot_pos <- -3
        }
        else if((gene_start >= TE_start && (gene_start-400) <= TE_start) || (TE_end <= gene_start && TE_end >= (gene_start-400))){
          plot_pos <- -4
        }
        else if((gene_start >= TE_start && (gene_start-500) <= TE_start) || (TE_end <= gene_start && TE_end >= (gene_start-500))) {
          plot_pos <- -5
        }
        else if((gene_start >= TE_start && (gene_start-600) <= TE_start) || (TE_end <= gene_start && TE_end >= (gene_start-600))){
          plot_pos <- -6
        }
        else if((gene_start >= TE_start && (gene_start-700) <= TE_start) || (TE_end <= gene_start && TE_end >= (gene_start-700))) {
          plot_pos <- -7
        }
        else if((gene_start >= TE_start && (gene_start-800) <= TE_start) || (TE_end <= gene_start && TE_end >= (gene_start-800))){
          plot_pos <- -8
        }
        else if((gene_start >= TE_start && (gene_start-900) <= TE_start) || (TE_end <= gene_start && TE_end >= (gene_start-900))) {
          plot_pos <- -9
        }
        else if((gene_start >= TE_start && (gene_start-1000) <= TE_start) || (TE_end <= gene_start && TE_end >= (gene_start-1000))){
          plot_pos <- -10
        }
        else if((gene_start >= TE_start && (gene_start-1100) <= TE_start) || (TE_end <= gene_start && TE_end >= (gene_start-1100))) {
          plot_pos <- -11
        }
        else if((gene_start >= TE_start && (gene_start-1200) <= TE_start) || (TE_end <= gene_start && TE_end >= (gene_start-1200))){
          plot_pos <- -12
        }
        else if((gene_start >= TE_start && (gene_start-1300) <= TE_start) || (TE_end <= gene_start && TE_end >= (gene_start-1300))) {
          plot_pos <- -13
        }
        else if((gene_start >= TE_start && (gene_start-1400) <= TE_start) || (TE_end <= gene_start && TE_end >= (gene_start-1400))){
          plot_pos <- -14
        }
        else if((gene_start >= TE_start && (gene_start-1500) <= TE_start) || (TE_end <= gene_start && TE_end >= (gene_start-1500))) {
          plot_pos <- -15
        }
        else if((gene_start >= TE_start && (gene_start-1600) <= TE_start) || (TE_end <= gene_start && TE_end >= (gene_start-1600))){
          plot_pos <- -16
        }
        else if((gene_start >= TE_start && (gene_start-1700) <= TE_start) || (TE_end <= gene_start && TE_end >= (gene_start-1700))) {
          plot_pos <- -17
        }
        else if((gene_start >= TE_start && (gene_start-1800) <= TE_start) || (TE_end <= gene_start && TE_end >= (gene_start-1800))){
          plot_pos <- -18
        }
        else if((gene_start >= TE_start && (gene_start-1900) <= TE_start) || (TE_end <= gene_start && TE_end >= (gene_start-1900))) {
          plot_pos <- -19
        }
        else if((gene_start >= TE_start && (gene_start-2000) <= TE_start) || (TE_end <= gene_start && TE_end >= (gene_start-2000))){
          plot_pos <- -20
        }
        else if((gene_start >= TE_start && (gene_start-2100) <= TE_start) || (TE_end <= gene_start && TE_end >= (gene_start-2100))) {
          plot_pos <- -21
        }
        else if((gene_start >= TE_start && (gene_start-2200) <= TE_start) || (TE_end <= gene_start && TE_end >= (gene_start-2200))){
          plot_pos <- -22
        }
        else if((gene_start >= TE_start && (gene_start-2300) <= TE_start) || (TE_end <= gene_start && TE_end >= (gene_start-2300))) {
          plot_pos <- -23
        }
        else if((gene_start >= TE_start && (gene_start-2400) <= TE_start) || (TE_end <= gene_start && TE_end >= (gene_start-2400))){
          plot_pos <- -24
        }
        else if((gene_start >= TE_start && (gene_start-2500) <= TE_start) || (TE_end <= gene_start && TE_end >= (gene_start-2500))) {
          plot_pos <- -25
        }
        else if((gene_start >= TE_start && (gene_start-2600) <= TE_start) || (TE_end <= gene_start && TE_end >= (gene_start-2600))){
          plot_pos <- -26
        }
        else if((gene_start >= TE_start && (gene_start-2700) <= TE_start) || (TE_end <= gene_start && TE_end >= (gene_start-2700))) {
          plot_pos <- -27
        }
        else if((gene_start >= TE_start && (gene_start-2800) <= TE_start) || (TE_end <= gene_start && TE_end >= (gene_start-2800))){
          plot_pos <- -28
        }
        else if((gene_start >= TE_start && (gene_start-2900) <= TE_start) || (TE_end <= gene_start && TE_end >= (gene_start-2900))) {
          plot_pos <- -29
        }
        else if((gene_start >= TE_start && (gene_start-3000) <= TE_start) || (TE_end <= gene_start && TE_end >= (gene_start-3000))){
          plot_pos <- -30
        }
        else if((gene_start >= TE_start && (gene_start-3100) <= TE_start) || (TE_end <= gene_start && TE_end >= (gene_start-3100))) {
          plot_pos <- -31
        }
        else if((gene_start >= TE_start && (gene_start-3200) <= TE_start) || (TE_end <= gene_start && TE_end >= (gene_start-3200))){
          plot_pos <- -32
        }
        else if((gene_start >= TE_start && (gene_start-3300) <= TE_start) || (TE_end <= gene_start && TE_end >= (gene_start-3300))) {
          plot_pos <- -33
        }
        else if((gene_start >= TE_start && (gene_start-3400) <= TE_start) || (TE_end <= gene_start && TE_end >= (gene_start-3400))){
          plot_pos <- -34
        }
        else if((gene_start >= TE_start && (gene_start-3500) <= TE_start) || (TE_end <= gene_start && TE_end >= (gene_start-3500))) {
          plot_pos <- -35
        }
        else if((gene_start >= TE_start && (gene_start-3600) <= TE_start) || (TE_end <= gene_start && TE_end >= (gene_start-3600))){
          plot_pos <- -36
        }
        else if((gene_start >= TE_start && (gene_start-3700) <= TE_start) || (TE_end <= gene_start && TE_end >= (gene_start-3700))) {
          plot_pos <- -37
        }
        else if((gene_start >= TE_start && (gene_start-3800) <= TE_start) || (TE_end <= gene_start && TE_end >= (gene_start-3800))){
          plot_pos <- -38
        }
        else if((gene_start >= TE_start && (gene_start-3900) <= TE_start) || (TE_end <= gene_start && TE_end >= (gene_start-3900))) {
          plot_pos <- -39
        }
        else if((gene_start >= TE_start && (gene_start-4000) <= TE_start) || (TE_end <= gene_start && TE_end >= (gene_start-4000))){
          plot_pos <- -40
        }
        else if((gene_start >= TE_start && (gene_start-4100) <= TE_start) || (TE_end <= gene_start && TE_end >= (gene_start-4100))) {
          plot_pos <- -41
        }
        else if((gene_start >= TE_start && (gene_start-4200) <= TE_start) || (TE_end <= gene_start && TE_end >= (gene_start-4200))){
          plot_pos <- -42
        }
        else if((gene_start >= TE_start && (gene_start-4300) <= TE_start) || (TE_end <= gene_start && TE_end >= (gene_start-4300))) {
          plot_pos <- -43
        }
        else if((gene_start >= TE_start && (gene_start-4400) <= TE_start) || (TE_end <= gene_start && TE_end >= (gene_start-4400))){
          plot_pos <- -44
        }
        else if((gene_start >= TE_start && (gene_start-4500) <= TE_start) || (TE_end <= gene_start && TE_end >= (gene_start-4500))) {
          plot_pos <- -45
        }
        else if((gene_start >= TE_start && (gene_start-4600) <= TE_start) || (TE_end <= gene_start && TE_end >= (gene_start-4600))){
          plot_pos <- -46
        }
        else if((gene_start >= TE_start && (gene_start-4700) <= TE_start) || (TE_end <= gene_start && TE_end >= (gene_start-4700))) {
          plot_pos <- -47
        }
        else if((gene_start >= TE_start && (gene_start-4800) <= TE_start) || (TE_end <= gene_start && TE_end >= (gene_start-4800))){
          plot_pos <- -48
        }
        else if((gene_start >= TE_start && (gene_start-4900) <= TE_start) || (TE_end <= gene_start && TE_end >= (gene_start-4900))) {
          plot_pos <- -49
        }
        else if((gene_start >= TE_start && (gene_start-5000) <= TE_start) || (TE_end <= gene_start && TE_end >= (gene_start-5000))){
          plot_pos <- -50
        }
      }
      # save the position to be plotted
      if(plot_pos != FALSE){
        data2plot <- c(data2plot,plot_pos)
        # save info about the TE in gene area
        a <- a+1
        TEinGENE[a,] <- c(TEsInRef[j,], as.character(DET_diapausa$Gene_ID[i]))
      }
    }
  }
}


## Quick visualization
summary(data2plot)
length(data2plot)
plot(table(data2plot))
hist(data2plot)
hist(data2plot,  breaks=500, xlim = c(-60,60), ylim=c(0,120))
plot(density(data2plot),xlim = c(-60,60), ylim=c(0.009,0.011)) 
