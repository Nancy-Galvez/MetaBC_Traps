#**Defining lineages at different clustering levels for each of the orders**
#**This script gets all scripts of the Species delimitation using 8 arthropods. For this, each of the ASV filtered dataset was used to generate an UPGMA tree with corrected distances under a F84 model, and based on this tree, all haplotypes were nested into clusterin levels (CLs) following the genetic similarity at different thresholds (0.5%, 1.5%, 3%, 5% and 7.5%), plus an additional threshold corresponding to the result of a species delimitation analyses conducted with the generalized mixed Yule-coalescent (GMYC).
#The analyses were performed using vegan (Oksanen et al., 2013), cluster, PMCMR, hier.part, ecodist, and betapart (Baselga & Orme, 2012).
####We used 8 groups: Arachnda, Coleoptera, Collembola, Diptera, Hemiptera, Hymenoptera, Lepidoptera, and Myriapoda. 

#**We get Part 1. Analysis of the UPGMAtree and GMYC each lineages**
library(igraph)
library(ape)
library(phangorn)
library(permute)
library(geiger)
library(lattice)
library(stats)
library(vegan)
library(base)
library(MASS)
library(paran)
library(gtools)
library(splits)
library(seqinr)
library(spam)
library(dotCall64)
library(grid)
library(maps)
library(fields)

##Source to scripts Diptera and data
source("Arachnida/SpeciesDelimitation/1UPGMAtree_AND_GMYC_each_lineage_ArachnidaStep1.R")
source("Coleoptera/SpeciesDelimitation/1UPGMAtree_AND_GMYC_each_lineage_ColeopteraStep1.R")
source("Collembola/SpeciesDelimitation/1UPGMAtree_AND_GMYC_each_lineage_CollembolaStep1.R")
source("Diptera/SpeciesDelimitation/1UPGMAtree_AND_GMYC_each_lineage_DipteraStep1.R")
source("Hemiptera/SpeciesDelimitation/1UPGMAtree_AND_GMYC_each_lineage_HemipteraStep1.R")
source("Hymenoptera/SpeciesDelimitation/1UPGMAtree_AND_GMYC_each_lineage_HymenopteraStep1.R")
source("Lepidoptera/SpeciesDelimitation/1UPGMAtree_AND_GMYC_each_lineage_LepidopteraStep1.R")                    
source("Myriapoda/SpeciesDelimitation/1UPGMAtree_AND_GMYC_each_lineage_MyriapodaStep1.R")


#**We get Part 2. Analysis to apply NODE.MIN get trees multiple each lineages**
#library(ape)
#library(phangorn)
#library(seqinr)
#library(stringr)

#source("Arachnida/to apply NODE.MIN_Get_trees_multiple_each_lineage_Arachnida.R") 
#source("Diptera/SpeciesDelimitation/to apply NODE.MIN_Get_trees_multiple_each_lineage_Diptera.R") 
#source("Collembola/SpeciesDelimitation/to apply NODE.MIN_Get_trees_multiple_each_lineage_Collembola.R") 
#source("Coleoptera/SpeciesDelimitation/to apply NODE.MIN_Get_trees_multiple_each_lineage_Coleoptera.R") 
#source("Hymenoptera/SpeciesDelimitation/to apply NODE.MIN_Get_trees_multiple_each_lineage_Hymenoptera.R") 
#source("Hemiptera/SpeciesDelimitation/to apply NODE.MIN_Get_trees_multiple_each_lineage_Hemiptera.R") 
#source("Lepidoptera/SpeciesDelimitation/to apply NODE.MIN_Get_trees_multiple_each_lineage_Lepidoptera.R") 
#source("Myriapoda/SpeciesDelimitation/to apply NODE.MIN_Get_trees_multiple_each_lineage_Myriapoda.R") 



#**END**
