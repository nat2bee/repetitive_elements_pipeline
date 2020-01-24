"""
## Date Created: 26-Nov-19
## Author: N. S. Araujo
## Script to test whether certain repetitive elements categories/families are enriched in a certain dataset in comparison to the entire genome.

E.g. 
Probability of getting TEs of the group DNA/TcMar-Tc1 around genes DE during diapause:
 >> n around gene= 519 / in n of all TEs around gene = length(TEinGENE$repeat_class.family) = 7755
 >> TEs of the group DNA/TcMar-Tc1 that we have in the entire genome = 50732 / in n of all TEs in the genome = length(all_TEs$repeat_class.family) = 818382
matrix with these data:
 >> matrix(c(519, 7755, 50732, 818382), 2, 2)
Statistics (one tailled):
 >> fisher.test(matrix(c(519, 7755, 50732, 818382), 2, 2), alternative='greater')
"""


## Data of interest for the stats
all_TEs_freq <- table(all_TEs$repeat_class.family)
allgenes_TEs_freq <- table(allgenes_TEinGENE$repeat_class.family)
DET_TEs_freq <- table(DET_TEinGENE$repeat_class.family)

all_TEs_size <- length(all_TEs$repeat_class.family)
allgenes_TEs_size <- length(allgenes_TEinGENE$repeat_class.family)
DET_TEs_size <- length(DET_TEinGENE$repeat_class.family)


## Loop for testing all repetitive elements categories in the genome library

for(i in 1:length(all_TEs_freq)){
  t_up <- 1
  t_down <- 1
  # if the TE family is also among the DET area
  if(names(all_TEs_freq[i]) %in% names(allgenes_TEs_freq)){
    TE <- names(all_TEs_freq[i])
    all_n <- as.numeric(all_TEs_freq[i])
    sample_n <- as.numeric(allgenes_TEs_freq[names(all_TEs_freq[i])])
  }
  # if it is not in the DET area
  else{
    TE <- names(all_TEs_freq[i])
    all_n <- as.numeric(all_TEs_freq[i])
    sample_n <- 0
  }
  # perform the statistical test 
  print(paste0("TE class '", TE, "'compared to the total genome (",sample_n,"vs",all_n,")"))
  t_up <- fisher.test(matrix(c(sample_n, allgenes_TEs_size, all_n, all_TEs_size), 2, 2), alternative='greater')
  t_down <- fisher.test(matrix(c(sample_n, allgenes_TEs_size, all_n, all_TEs_size), 2, 2), alternative='less')
  if(t_up$p.value <= 0.01){
    print(t_up)
  }
  if(t_down$p.value <= 0.01){
    print(t_down)
  }
}


for(i in 1:length(allgenes_TEs_freq)){
  t_up <- 1
  t_down <- 1
  # if the TE family is also among the DET area
  if(names(allgenes_TEs_freq[i]) %in% names(DET_TEs_freq)){
    TE <- names(allgenes_TEs_freq[i])
    all_n <- as.numeric(allgenes_TEs_freq[i])
    sample_n <- as.numeric(DET_TEs_freq[names(allgenes_TEs_freq[i])])
  }
  # if it is not in the DET area
  else{
    TE <- names(allgenes_TEs_freq[i])
    all_n <- as.numeric(allgenes_TEs_freq[i])
    sample_n <- 0
  }
  # perform the statistical test 
  print(paste0("TE class '", TE, "'compared to the total genome (",sample_n,"vs",all_n,")"))
  t_up <- fisher.test(matrix(c(sample_n, DET_TEs_size, all_n, allgenes_TEs_size), 2, 2), alternative='greater')
  t_down <- fisher.test(matrix(c(sample_n, DET_TEs_size, all_n, allgenes_TEs_size), 2, 2), alternative='less')
  if(t_up$p.value <= 0.01){
    print(t_up)
  }
  if(t_down$p.value <= 0.01){
    print(t_down)
  }
}
