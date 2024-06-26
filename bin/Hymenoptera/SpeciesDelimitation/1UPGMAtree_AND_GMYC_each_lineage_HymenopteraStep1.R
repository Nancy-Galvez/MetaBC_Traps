#**Defining Hymenoptera lineage at different clustering levels: STEP 1** 

#**This script gets Specie delimitation using **Hymenoptera order**. For this, the ASV filtered dataset was used to generate an UPGMA tree with corrected distances under a F84 model, and based on this tree, all haplotypes were nested into clusterin levels (CLs) following the genetic similarity at different thresholds (0.5%, 1.5%, 3%, 5% and 7.5%), plus an additional threshold corresponding to the result of a specie delimitation analyses conducted with the generalized mixed Yule-coalescent (GMYC).

#PART 1. Distance matrix and UPGMA 
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


Hymenoptera <- read.FASTA("../../../genetic/Orders9&Etiquetas/Hymenoptera_161_33SOPAS_UNIQUES_ExpeRelabel420_20_26_1_alignment.fasta")  #leer fasta
dist.dna(Hymenoptera, model = "F84", as.matrix = T)->dist_Hymenoptera  #calcular matriz distancias por pares
Hymenoptera_UPGMA <- upgma(dist_Hymenoptera)  #hacer árbol
plot(Hymenoptera_UPGMA)

#PART 2: GMYC

#Those funsions are in bin
source("../../gmyc.pkg.0.9.6.R")
source("../../Powell_supplemental_script.R")

#**single gmyc threshold** 
Hymenoptera_UPGMA.GMYC.singleThreshold<-gmyc.edit(Hymenoptera_UPGMA, method="s") #hacer GMYC
#resumen de los resultados
summary(Hymenoptera_UPGMA.GMYC.singleThreshold) #Mostrar resumen GMYC

sink("../../Hymenoptera/SpeciesDelimitation/Hymenoptera_UPGMA.GMYC.singleThreshold.log", type=c("output","message"), split=TRUE) #guardar resultado del GMYC (resumen)
#MIT.test.s o MIT.test.m para ver todos los datos
summary(Hymenoptera_UPGMA.GMYC.singleThreshold)
sink()

#We used the result of the threshold time: e.g. "0.0292959" in the next analisys.

#**Following step 2**




