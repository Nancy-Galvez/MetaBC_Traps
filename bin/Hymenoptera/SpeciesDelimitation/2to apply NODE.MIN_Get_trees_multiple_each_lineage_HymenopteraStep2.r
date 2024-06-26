#**Defining Hemiptera lineage at different clustering levels: STEP 2 and 3** 

#**This script gets Specie delimitation using **Hemiptera order**. For this, the ASV filtered dataset was used to generate an UPGMA tree with corrected distances under a F84 model, and based on this tree, all haplotypes were nested into clusterin levels (CLs) following the genetic similarity at different thresholds (0.5%, 1.5%, 3%, 5% and 7.5%), plus an additional threshold corresponding to the result of a specie delimitation analyses conducted with the generalized mixed Yule-coalescent (GMYC).

#**PART 1. Distance Matrix and UPGMA**
library(ape)
library(phangorn)
library(seqinr)
library(stringr)

Hemiptera <- read.FASTA("../../../genetic/Orders9&Etiquetas/Hymenoptera_161_33SOPAS_UNIQUES_ExpeRelabel420_20_26_1_alignment.fasta")  #leer fasta
dist.dna(Hemiptera, model = "F84", as.matrix = T)->dist_Hemiptera  #calcular matriz distancias por pares
Hemiptera_UPGMA <- upgma(dist_Hemiptera)  #hacer arbol
plot(Hemiptera_UPGMA)

#PART 3: To obtain subtrees at multihierarchical levels.
#This funcion are in bin. NOTA Nos muestra un warning, aunque sigue estando bien.

source("../../NODE.MIN_Get_trees_multiple.R")

#**Hemiptera**
#We used the result of the threshold time: e.g. "0.0292959" in the next analisys of each lineages.
read.table ("../../../genetic/Orders9&Etiquetas/Hymenoptera_161_33SOPAS_UNIQUES_ExpeRelabel420_20_26_1_alignment.fasta", header=FALSE)->Hemiptera.fasta
#read.tree("Hemiptera_All_s2_sinSTOPS_revised.newick")->Hemiptera_SHORT2.tree
GMYC.thershold<-0.01095415
		
NODE.MIN(Hemiptera.fasta,Hemiptera_UPGMA,0.005, 0.01, 0.015, 0.02, 0.025, 0.03, 0.04, 0.05, 0.06, 0.075, 0.075, 0.075, GMYC.thershold, print.subtrees = "NO", print.subtrees.fasta = "NO")

#**END**

