#**INITIAL STEPS**

#**Community diversity and composition at the 3% Clustering level of the Coleoptera order**

#####In excel remove the simbol ## from the names of the table and rename the samples acccording to the code used in the gradient e.g. GRA_S10_D_F_A10
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

#**TABLES AND COMMUNITY MATRIXES**
#####open table with names including Region and habitat parameters
s2_raw_all <- read.table("../../genetic/Data_in/Coleoptera/s2_raw_all_Coleoptera_threshold.txt", header=TRUE)
dim(s2_raw_all)

#####remove additional columns and leave only names (of haplotipes), samples and taxa (and threshold in this case)
s2_raw_all[,c(1:44,46)]->s2_raw
dim(s2_raw) ##49 samples = 48 plus 1 neg (the second neg from DOM_REPS is not there because all 0)
colnames(s2_raw)

#####Applying the conservative threshold (this is a binary column)
s2_raw[which(s2_raw$conservative_threshold == "1"),]->s2_raw_threshold 
s2_raw_threshold [,1:44]->s2_raw_threshold ##remove threshold col
dim(s2_raw_threshold)
colnames(s2_raw_threshold)

#####loop to create the new matrix combining haplotype by otu pertenency, i.e. submatrix by limit
#**Coleoptera**
unique (s2_raw_threshold$limite0.03)->levels_limite0.03

data.frame()->s2_raw_Coleoptera_limite0.03

for (i in 1:length(unique (s2_raw_threshold$limite0.03)))
{
  levels_limite0.03[i]->level0.03
  s2_raw_threshold[which(s2_raw_threshold$limite0.03==level0.03),]->subcom_level_names0.03
  subcom_level_names0.03[,c(2:30)]->subcom_level0.03  #delete names, level and also the negative column
  colSums(subcom_level0.03)->sum0.03
  as.data.frame(sum0.03)->sum0.03
  t(sum0.03)->sum0.03
  row.names(sum0.03)<-subcom_level_names0.03[1,1] #keep the name of the first haplotype
  rbind(s2_raw_Coleoptera_limite0.03,sum0.03)->s2_raw_Coleoptera_limite0.03
}

##**transform in present/absence table**
s2_raw_Coleoptera_limite0.03->s2_raw_Coleoptera_limite0.03
s2_raw_Coleoptera_limite0.03[s2_raw_Coleoptera_limite0.03>1]<-1 ##transform in present/absence table 

##**checking if there is any row with no presence**
s2_raw_Coleoptera_limite0.03[,1:29]->data0.03
rowSums(data0.03)
length(which(rowSums(data0.03)!=0))
length(which(rowSums(data0.03)==0))

##**Community matrixes (samples in rows and h in cols).**
##**Coleoptera**
t(s2_raw_Coleoptera_limite0.03)->t_s2_f4_Coleoptera_limite0.03 ##trasp
t_s2_f4_Coleoptera_limite0.03[1:29,]->community_Coleoptera_limite0.03 #NOTA_Nancy: Este numero es importante. Colocar exactamente el numero de "s2_f4[,2:52]->data0.03".
colnames(t_s2_f4_Coleoptera_limite0.03)<-community_Coleoptera_limite0.03[1,] 
as.data.frame(community_Coleoptera_limite0.03)->community_Coleoptera0.03 ##trasp including col and row names
####community_Acari[-49,]->community_Coleoptera ##removing neg
dim(community_Coleoptera0.03)
community_Coleoptera0.03[order(row.names(community_Coleoptera0.03)),]->community_Coleoptera0.03 ##order samples
write.table (community_Coleoptera0.03, file="../../genetic/Data_out/Coleoptera/Coleoptera3P/community_Coleoptera0.03.txt") ##this is necessary for the format, not able to solve in other way
read.table ("../../genetic/Data_out/Coleoptera/Coleoptera3P/community_Coleoptera0.03.txt")->community_Coleoptera0.03


##**submatrixes by Traps in Nevado Toluca**
dim(community_Coleoptera0.03)
community_Coleoptera0.03[which(str_extract (row.names(community_Coleoptera0.03), "_T_") %in% "_T_"),]->community_Coleoptera_Traps0.03
dim(community_Coleoptera_Traps0.03)
community_Coleoptera_Traps0.03[,which(colSums(community_Coleoptera_Traps0.03)!=0)]->community_Coleoptera_Traps0.03 ##to remove no data colums
dim(community_Coleoptera_Traps0.03)
write.table (community_Coleoptera_Traps0.03, file="../../genetic/Data_out/Coleoptera/Coleoptera3P/community_Coleoptera_Traps0.03.txt") ##this is necessary for the format, not able to solve in other way
read.table ("../../genetic/Data_out/Coleoptera/Coleoptera3P/community_Coleoptera_Traps0.03.txt")->community_Coleoptera_Traps0.03


##**Generating a general table with names and habitat parameters.**
##BY Traps
##**Generating a general table with names and habitat parameters**
row.names(community_Coleoptera_Traps0.03)->sample_names_Mountain1_0.03
as.data.frame(sample_names_Mountain1_0.03)->sample_names_Mountain1_0.03
sample_names_Mountain1_0.03 %>% separate(sample_names_Mountain1_0.03, c("Conservation1","Mountain1","Traps","ID"), sep="_",remove=FALSE)->general_sample_Mountain1Traps0.03
general_sample_Mountain1Traps0.03
general_sample_Mountain1Traps0.03 %>% unite(Mountain1andTraps, Mountain1,Traps, sep="_",remove=FALSE)->general_sample_Mountain1Traps0.03 ####generating a variable combining layer and habitat
general_sample_Mountain1Traps0.03
write.table(general_sample_Mountain1Traps0.03, file="../../genetic/Data_out/Coleoptera/Coleoptera3P/general_sample_Mountain1Traps0.03.txt") ##this is the only way I found to be able to work later
read.table("../../genetic/Data_out/Coleoptera/Coleoptera3P/general_sample_Mountain1Traps0.03.txt",header=TRUE)->general_sample_Mountain1Traps0.03



#**HAPLOTYPE RICHNESS TABLES, PLOTS AND ANALYSES by TrapsS**
##**Coleoptera**
as.matrix(community_Coleoptera_Traps0.03)->community_Coleoptera_Traps0.03
row.names(community_Coleoptera_Traps0.03)->sample_names_Traps0.03
dim(community_Coleoptera_Traps0.03)->dims_Traps0.03
dims_Traps0.03
.rowSums (community_Coleoptera_Traps0.03,dims_Traps0.03[1],dims_Traps0.03[2])->sample_richness_Traps0.03 ##summatory by rows
rbind(sample_names_Traps0.03,sample_richness_Traps0.03)->richness_Traps0.03
t(richness_Traps0.03)->richness_Traps0.03
colnames(richness_Traps0.03)<-c("sample_names_Traps0.03","sample_richness_Traps0.03")
richness_Traps0.03
as.data.frame(richness_Traps0.03)->richness_Traps0.03

##**Generating variables with Traps. En la tabla rishnees. Generar variable de montana y sitio.** 
richness_Traps0.03 %>% separate(sample_names_Traps0.03, c("Conservation1","Mountain1","Traps","ID"), sep="_",remove=FALSE)->richness_Traps0.03
richness_Traps0.03
richness_Traps0.03 %>% unite(Mountain1Traps, Mountain1, Traps, sep="_",remove=FALSE)->richness_Traps0.03 ##generating a variable combining layer and habitat
richness_Traps0.03

##**Generating variables with in total Traps_C.**
#richness_Traps0.03 %>% separate(sample_names_Traps0.03, c("Conservation","Mountain1","Traps","ID"), sep="_",remove=FALSE)->richness_TrapsC
#richness_TrapsC
#richness_TrapsC %>% unite(ConservationMountain1, Conservation, Mountain1, sep="_",remove=FALSE)->richness_TrapsC ##generating a variable combining layer and habitat
#richness_TrapsC

##**BY Traps**
write.table(richness_Traps0.03, file="../../genetic/Data_out/Coleoptera/Coleoptera3P/richness_Traps0.03_Coleoptera.txt") ##this is the only way I found to be able to work later
read.table("../../genetic/Data_out/Coleoptera/Coleoptera3P/richness_Traps0.03_Coleoptera.txt",header=TRUE)->richness_Traps0.03

##**BY Traps_C general**
#write.table(richness_TrapsC, file="../../genetic/Data_out/Coleoptera/Coleoptera3P/richness_TrapsC_Coleoptera.txt") ##this is the only way I found to be able to work later
#read.table("../../genetic/Data_out/Coleoptera/Coleoptera3P/richness_TrapsC_Coleoptera.txt",header=TRUE)->richness_TrapsC

##**General plot of richness by sample in Traps**
barplot(richness_Traps0.03$sample_richness_Traps0.03,col=richness_Traps0.03$Mountain1Traps,names.arg= richness_Traps0.03$sample_names_Traps0.03,las=2,cex.names=0.5, ylab="richness_Traps0.03", main="H richness_Traps Coleoptera_0.03")
richness_Traps0.03 %>% group_by(Mountain1Traps) %>% summarise(mean(sample_richness_Traps0.03))

##**min, max, ds Summarise**
richness_Traps0.03 %>% group_by(Mountain1Traps) %>% summarise(min(sample_richness_Traps0.03))
richness_Traps0.03 %>% group_by(Mountain1Traps) %>% summarise(max(sample_richness_Traps0.03))
richness_Traps0.03 %>% group_by(Mountain1Traps) %>% summarise(sd(sample_richness_Traps0.03))

##**General mean, min, max, ds by sample in Traps_richness_TrapsC**
#richness_TrapsC %>% group_by(ConservationMountain1) %>% summarise(mean(sample_richness_Traps0.03))
#richness_TrapsC %>% group_by(ConservationMountain1) %>% summarise(min(sample_richness_Traps0.03))
#richness_TrapsC %>% group_by(ConservationMountain1) %>% summarise(max(sample_richness_Traps0.03))
#richness_TrapsC %>% group_by(ConservationMountain1) %>% summarise(sd(sample_richness_Traps0.03))

##**Global richness by Traps.**  
plot(richness_Traps0.03$Mountain1Traps,richness_Traps0.03$sample_richness_Traps0.03,ylab="richness_Traps0.03", ylim=c(0,6), cex=1.4, cex.axis=2.3, lwd=2.5)
kruskal.test(sample_richness_Traps0.03 ~ Mountain1Traps, data = richness_Traps0.03)
#posthoc.kruskal.nemenyi.test(x=richness_Traps0.03$sample_richness_Traps0.03, g=richness_Traps0.03$Mountain1Traps, method="Bonferroni")
PMCMRplus::kwAllPairsConoverTest(x=richness_Traps0.03$sample_richness_Traps0.03, g=richness_Traps0.03$Mountain1Traps, p.adjust.method = "single-step")
##Comparison of each group against. 
#text(x=c(1,2,3,4.1), y=(49), labels=c("ab","a","ab","b"), cex=1.4)
#text(x=c(1,2,3,4), y=(2), labels=c("ab","a","ab","b"), cex=1.5)
text(x=4.3, y=0.5, labels="ns", cex=2)
mtext(c("lineages 3%"), side = 3, col = "black", line = 1, cex = 2)



#**BETADIVERSITY ORDINATIONS by Trapss** 
#**Coleoptera**
##beta general

beta.multi(community_Coleoptera_Traps0.03, index.family="sorensen")

##turnover by pairs, nmds, anosim
beta.pair(community_Coleoptera_Traps0.03, index.family="sorensen")->beta.pair  ##betadiversity by pair of communities using sorensen on the precense/absence data, with estimation of turnover and nestedness datamatrixes simultaneously
metaMDS (beta.pair$beta.sim)->MDSbetasim_T0.03 ##NMDS

plot (MDSbetasim_T0.03, main="Coleoptera_Traps_0.03") 
x<- MDSbetasim_T0.03$points[,1]
y<- MDSbetasim_T0.03$points[,2]
text(x, y, pos = 1, cex=0.7, labels = row.names (community_Coleoptera_Traps0.03))

plot (MDSbetasim_T0.03, main="Coleoptera_Traps_0.03")
with(general_sample_Mountain1Traps0.03,ordispider(MDSbetasim_T0.03, Traps, label=T, col="blue"))

plot (MDSbetasim_T0.03, xlim=c(-0.4, 0.4), ylim=c(-0.4, 0.4), cex.axis=1.4, cex=1.2, cex.lab=1.4, main="Coleoptera_Traps_0.03")
with(general_sample_Mountain1Traps0.03,ordispider(MDSbetasim_T0.03, Traps, label=T, cex.lab=0.9, col= c("#153a7b", "#eaa22f", "#97518b", "#81b9d0")))

##**Anosim**
anosim(beta.pair$beta.sim, general_sample_Mountain1Traps0.03$Traps, permutations=999)

plot (MDSbetasim_T0.03, xlim=c(-0.4, 0.4), ylim=c(-0.4, 0.4), cex=2, cex.lab=1, cex.axis=2.3, lwd=4.5)
with(general_sample_Mountain1Traps0.03,ordispider(MDSbetasim_T0.03, Traps, label=T, cex.lab=1, col= c("#153a7b", "#eaa22f", "#97518b", "#81b9d0"), lwd=4.5))
##I put letter "r" in cursive and r2 value
mylabel = bquote(italic(r)^2)
text(x=0.20, y=-0.35, labels = mylabel, cex=2)
text(x=0.42, y=-0.35, labels="=0.28 *", cex=2)
mtext(c("lineages 3%"), side = 3, col = "black", line = 1, cex = 2)


#**END TRAPS**


###################################
##**submatrixes by Groups in Nevado Toluca**
dim(community_Coleoptera0.03)
community_Coleoptera0.03[which(str_extract (row.names(community_Coleoptera0.03), "G_") %in% "G_"),]->community_Coleoptera_Groups0.03
dim(community_Coleoptera_Groups0.03)
community_Coleoptera_Groups0.03[,which(colSums(community_Coleoptera_Groups0.03)!=0)]->community_Coleoptera_Groups0.03 ##to remove no data colums
dim(community_Coleoptera_Groups0.03)
write.table (community_Coleoptera_Groups0.03, file="../../genetic/Data_out/Coleoptera/Coleoptera3P/community_Coleoptera_Groups0.03.txt") ##this is necessary for the format, not able to solve in other way
read.table ("../../genetic/Data_out/Coleoptera/Coleoptera3P/community_Coleoptera_Groups0.03.txt")->community_Coleoptera_Groups0.03


##**Generating a general table with names and habitat parameters.**
##BY Groups
##**Generating a general table with names and habitat parameters**
row.names(community_Coleoptera_Groups0.03)->sample_names_Mountain2_0.03
as.data.frame(sample_names_Mountain2_0.03)->sample_names_Mountain2_0.03
sample_names_Mountain2_0.03 %>% separate(sample_names_Mountain2_0.03, c("Conservation2","Mountain2","Groups","ID"), sep="_",remove=FALSE)->general_sample_Mountain2Groups0.03
general_sample_Mountain2Groups0.03
general_sample_Mountain2Groups0.03 %>% unite(Mountain2andGroups, Mountain2,Groups, sep="_",remove=FALSE)->general_sample_Mountain2Groups0.03 ####generating a variable combining layer and habitat
general_sample_Mountain2Groups0.03
write.table(general_sample_Mountain2Groups0.03, file="../../genetic/Data_out/Coleoptera/Coleoptera3P/general_sample_Mountain2Groups0.03.txt") ##this is the only way I found to be able to work later
read.table("../../genetic/Data_out/Coleoptera/Coleoptera3P/general_sample_Mountain2Groups0.03.txt",header=TRUE)->general_sample_Mountain2Groups0.03



#**HAPLOTYPE RICHNESS TABLES, PLOTS AND ANALYSES by GroupsS**
##**Coleoptera**
as.matrix(community_Coleoptera_Groups0.03)->community_Coleoptera_Groups0.03
row.names(community_Coleoptera_Groups0.03)->sample_names_Groups0.03
dim(community_Coleoptera_Groups0.03)->dims_Groups0.03
dims_Groups0.03
.rowSums (community_Coleoptera_Groups0.03,dims_Groups0.03[1],dims_Groups0.03[2])->sample_richness_Groups0.03 ##summatory by rows
rbind(sample_names_Groups0.03,sample_richness_Groups0.03)->richness_Groups0.03
t(richness_Groups0.03)->richness_Groups0.03
colnames(richness_Groups0.03)<-c("sample_names_Groups0.03","sample_richness_Groups0.03")
richness_Groups0.03
as.data.frame(richness_Groups0.03)->richness_Groups0.03

##**Generating variables with Groups. En la tabla rishnees. Generar variable de montana y sitio.** 
richness_Groups0.03 %>% separate(sample_names_Groups0.03, c("Conservation2","Mountain2","Groups","ID"), sep="_",remove=FALSE)->richness_Groups0.03
richness_Groups0.03
richness_Groups0.03 %>% unite(Mountain2Groups, Mountain2, Groups, sep="_",remove=FALSE)->richness_Groups0.03 ##generating a variable combining layer and habitat
richness_Groups0.03

##**Generating variables with in total Groups_C.**
#richness_Groups0.03 %>% separate(sample_names_Groups0.03, c("Conservation","Mountain1","Groups","ID"), sep="_",remove=FALSE)->richness_GroupsC
#richness_GroupsC
#richness_GroupsC %>% unite(ConservationMountain1, Conservation, Mountain1, sep="_",remove=FALSE)->richness_GroupsC ##generating a variable combining layer and habitat
#richness_GroupsC

##**BY Groups**
write.table(richness_Groups0.03, file="../../genetic/Data_out/Coleoptera/Coleoptera3P/richness_Groups0.03_Coleoptera.txt") ##this is the only way I found to be able to work later
read.table("../../genetic/Data_out/Coleoptera/Coleoptera3P/richness_Groups0.03_Coleoptera.txt",header=TRUE)->richness_Groups0.03

##**BY Groups_C general**
#write.table(richness_GroupsC, file="../../genetic/Data_out/Coleoptera/Coleoptera3P/richness_GroupsC_Coleoptera.txt") ##this is the only way I found to be able to work later
#read.table("../../genetic/Data_out/Coleoptera/Coleoptera3P/richness_GroupsC_Coleoptera.txt",header=TRUE)->richness_GroupsC

##**General plot of richness by sample in Groups**
barplot(richness_Groups0.03$sample_richness_Groups0.03,col=richness_Groups0.03$Mountain2Groups,names.arg= richness_Groups0.03$sample_names_Groups0.03,las=2,cex.names=0.5, ylab="richness_Groups0.03", main="H richness_Groups Coleoptera_0.03")
richness_Groups0.03 %>% group_by(Mountain2Groups) %>% summarise(mean(sample_richness_Groups0.03))

##**min, max, ds Summarise**
richness_Groups0.03 %>% group_by(Mountain2Groups) %>% summarise(min(sample_richness_Groups0.03))
richness_Groups0.03 %>% group_by(Mountain2Groups) %>% summarise(max(sample_richness_Groups0.03))
richness_Groups0.03 %>% group_by(Mountain2Groups) %>% summarise(sd(sample_richness_Groups0.03))

##**General mean, min, max, ds by sample in Groups_richness_GroupsC**
#richness_GroupsC %>% group_by(ConservationMountain1) %>% summarise(mean(sample_richness_Groups0.03))
#richness_GroupsC %>% group_by(ConservationMountain1) %>% summarise(min(sample_richness_Groups0.03))
#richness_GroupsC %>% group_by(ConservationMountain1) %>% summarise(max(sample_richness_Groups0.03))
#richness_GroupsC %>% group_by(ConservationMountain1) %>% summarise(sd(sample_richness_Groups0.03))


##**Global richness by Groups.**  
plot(richness_Groups0.03$Mountain2Groups,richness_Groups0.03$sample_richness_Groups0.03,ylab="richness_Groups_0.03_Coleoptera", cex=1.4, cex.axis=2.3, lwd=2.5, ylim=c(0,10))
kruskal.test(sample_richness_Groups0.03 ~ Mountain2Groups, data = richness_Groups0.03)
#PMCMRplus::kwAllPairsNemenyiTest(x=richness_Groups_h_Coleoptera$sample_richness_Groups_h, g=richness_Groups_h_Coleoptera$Mountain2Groups, p.adjust.method = "single-step")
PMCMRplus::kwAllPairsConoverTest(x=richness_Groups0.03$sample_richness_Groups0.03, g=richness_Groups0.03$Mountain2Groups, p.adjust.method = "single-step")
#PMCMRplus::kwAllPairsDunnTest(x=richness_Groups_h_Coleoptera$sample_richness_Groups_h, g=richness_Groups_h_Coleoptera$Mountain2Groups, p.adjust.method="bonferroni")
#kwAllPairsNemenyiTest(x=richness_Groups_h_Coleoptera$sample_richness_Groups_h, g=richness_Groups_h_Coleoptera$Mountain2Groups)
##Comparison of each group against. 
#text(x=c(1,2,3,4,5), y=1.4, labels=c("ab","ab","ab","a","b"), cex=1.5)
text(x=5.4, y=1, labels="ns", cex=2)
mtext(c("lineages 3%"), side = 3, col = "black", line = 1, cex = 2)

#**BETADIVERSITY ORDINATIONS by Groupss** 
#**Coleoptera**
##beta general

beta.multi(community_Coleoptera_Groups0.03, index.family="sorensen")

##turnover by pairs, nmds, anosim
beta.pair(community_Coleoptera_Groups0.03, index.family="sorensen")->beta.pair  ##betadiversity by pair of communities using sorensen on the precense/absence data, with estimation of turnover and nestedness datamatrixes simultaneously
metaMDS (beta.pair$beta.sim)->MDSbetasim_G0.03 ##NMDS

plot (MDSbetasim_G0.03, main="Coleoptera_Groups_0.03") 
x<- MDSbetasim_G0.03$points[,1]
y<- MDSbetasim_G0.03$points[,2]
text(x, y, pos = 1, cex=0.7, labels = row.names (community_Coleoptera_Groups0.03))

##**By SITES quitar**
##**repeating after removing outlayers from matrix and general habitat table**
community_Coleoptera_Groups0.03[-which(row.names(community_Coleoptera_Groups0.03) %in% c("G_B_Sa_55TRE137a", "G_E_Pab_55TRE137")),]->community_Coleoptera_Groups0.03_sinoutlayer
community_Coleoptera_Groups0.03_sinoutlayer[,which(colSums(community_Coleoptera_Groups0.03_sinoutlayer)!=0)]->community_Coleoptera_Groups0.03_sinoutlayer ##to remove no data colums
dim(community_Coleoptera_Groups0.03_sinoutlayer)  

general_sample_Mountain2Groups0.03[-which(general_sample_Mountain2Groups0.03$sample_names %in% c("G_B_Sa_55TRE137a", "G_E_Pab_55TRE137")),]->general_sample_Groups0.03_sinoutlayer

beta.pair(community_Coleoptera_Groups0.03_sinoutlayer, index.family="sorensen")->beta.pair  ##betadiversity by pair of communities using sorensen on the precense/absence data, with estimation of turnover and nestedness datamatrixes simultaneously
metaMDS (beta.pair$beta.sim)->MDSbetasim_G10.03 ##NMDS

plot (MDSbetasim_G10.03, main="Coleoptera_Groups0.03") 
x<- MDSbetasim_G10.03$points[,1]
y<- MDSbetasim_G10.03$points[,2]
text(x, y, pos = 1, cex=0.7, labels = row.names (community_Coleoptera_Groups0.03_sinoutlayer))


plot (MDSbetasim_G10.03, main="Coleoptera_Groups0.03")
with(general_sample_Groups0.03_sinoutlayer,ordispider(MDSbetasim_G10.03, Mountain2, label=T, col="blue"))

plot (MDSbetasim_G10.03, xlim=c(-0.5, 0.5), ylim=c(-0.4, 0.4), cex.axis=1.4, cex=1.2, cex.lab=1.4, main="Coleoptera_Groups0.03")
with(general_sample_Groups0.03_sinoutlayer,ordispider(MDSbetasim_G10.03, Mountain2, label=T, cex.lab=0.9, col= c("#153a7b", "#eaa22f", "#97518b", "#81b9d0", "#c1d6b1")))

##**Anosim**
anosim(beta.pair$beta.sim, general_sample_Groups0.03_sinoutlayer$Mountain2, permutations=999)

plot (MDSbetasim_G10.03, xlim=c(-0.4, 0.4), ylim=c(-0.4, 0.4), cex=2, cex.lab=1, cex.axis=2.3, lwd=4.5)
with(general_sample_Groups0.03_sinoutlayer,ordispider(MDSbetasim_G10.03, Mountain2, label=T, cex.lab=1, col= c("#153a7b", "#eaa22f", "#97518b", "#81b9d0", "#c1d6b1"), lwd=4.5))
##I put letter "r" in cursive and r2 value
mylabel = bquote(italic(r)^2)
text(x=0.20, y=-0.35, labels = mylabel, cex=2)
text(x=0.45, y=-0.35, labels="=0.19 ns", cex=2)
mtext(c("lineages 3%"), side = 3, col = "black", line = 1, cex = 2)

#**END**
