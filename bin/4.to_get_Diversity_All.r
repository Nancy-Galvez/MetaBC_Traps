#**Community diversity and composition at the haplotype, 3% and 5% lineages for each of the five taxonomic orders studied**

#####This script gets all scripts of diversity using 8 arthropods order at multi-hierarchical levels 
#####We used 6 groups: Arachnda, Coleoptera, Collembola, Diptera, Hemiptera, and Hymenoptera.

library(stats)
library(base)
library(dplyr)
library(tidyr)
library(knitr)
library(PMCMR)
library(rcdd) 
library(vegan)
library(betapart) 
library(stringr)
library(permute)
library(lattice)


#**We get the diversity analysis at haplotypes, 3% and 5% lineages in Diptera**
##Source to scripts Diptera and data
source("Diptera/Com_matrixes_data_exploration_Diptera_h_TG.R") 
source("Diptera/Com_matrixes_data_exploration_Diptera_3P_TG.R") 
source("Diptera/Com_matrixes_data_exploration_Diptera_5P_TG.R") 

#**We get the diversity analysis at haplotypes, 3% and 5% lineages in Collembola**
##Source to scripts Collembola and data
source("Collembola/Com_matrixes_data_exploration_Collembola_h_TG.R") 
source("Collembola/Com_matrixes_data_exploration_Collembola_3P_TG.R") 
source("Collembola/Com_matrixes_data_exploration_Collembola_5P_TG.R") 

#**We get the diversity analysis at haplotypes, 3% and 5% lineages in Arachnida**
##Source to scripts Arachnida and data
source("Arachnida/Com_matrixes_data_exploration_Arachnida_h_TG.R") 
source("Arachnida/Com_matrixes_data_exploration_Arachnida_3P_TG.R") 
source("Arachnida/Com_matrixes_data_exploration_Arachnida_5P_TG.R")

#**We get the diversity analysis at haplotypes, C3% and 5% lineages in Hemiptera**
##Source to scripts Hemiptera and data
source("Hemiptera/Com_matrixes_data_exploration_Hemiptera_h_TG.R") 
source("Hemiptera/Com_matrixes_data_exploration_Hemiptera_3P_TG.R") 
source("Hemiptera/Com_matrixes_data_exploration_Hemiptera_5P_TG.R") 

#**We get the diversity analysis at haplotypes, 3% and 5% lineages in Hymenoptera**
##Source to scripts Hymenoptera and data
source("Hymenoptera/Com_matrixes_data_exploration_Hymenoptera_h_TG.R") 
source("Hymenoptera/Com_matrixes_data_exploration_Hymenoptera_3P_TG.R") 
source("Hymenoptera/Com_matrixes_data_exploration_Hymenoptera_5P_TG.R") 

#**We get the diversity analysis at haplotypes, 3% and 5% lineages in Coleoptera**
##Source to scripts Coleoptera and data
source("Coleoptera/Com_matrixes_data_exploration_Coleoptera_h_TG.R") 
source("Coleoptera/Com_matrixes_data_exploration_Coleoptera_3P_TG.R") 
source("Coleoptera/Com_matrixes_data_exploration_Coleoptera_5P_TG.R") 

#**END**
