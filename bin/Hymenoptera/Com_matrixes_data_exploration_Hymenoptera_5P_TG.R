#**INITIAL STEPS**

#**Community diversity and composition at the 5% Clustering level of the Hymenoptera order**

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
s2_raw_all <- read.table("../../genetic/Data_in/Hymenoptera/s2_raw_all_Hymenoptera_threshold.txt", header=TRUE)
dim(s2_raw_all)

#####remove additional columns and leave only names (of haplotipes), samples and taxa (and threshold in this case)
s2_raw_all[,c(1:42,44)]->s2_raw
dim(s2_raw) ##49 samples = 48 plus 1 neg (the second neg from DOM_REPS is not there because all 0)
colnames(s2_raw)

#####Applying the conservative threshold (this is a binary column)
s2_raw[which(s2_raw$conservative_threshold == "1"),]->s2_raw_threshold 
s2_raw_threshold [,1:42]->s2_raw_threshold ##remove threshold col
dim(s2_raw_threshold)
colnames(s2_raw_threshold)

#####loop to create the new matrix combining haplotype by otu pertenency, i.e. submatrix by limit
#**Hymenoptera**
unique (s2_raw_threshold$limite0.05)->levels_limite0.05

data.frame()->s2_raw_Hymenoptera_limite0.05

for (i in 1:length(unique (s2_raw_threshold$limite0.05)))
{
  levels_limite0.05[i]->level0.05
  s2_raw_threshold[which(s2_raw_threshold$limite0.05==level0.05),]->subcom_level_names0.05
  subcom_level_names0.05[,c(2:28)]->subcom_level0.05  #delete names, level and also the negative column
  colSums(subcom_level0.05)->sum0.05
  as.data.frame(sum0.05)->sum0.05
  t(sum0.05)->sum0.05
  row.names(sum0.05)<-subcom_level_names0.05[1,1] #keep the name of the first haplotype
  rbind(s2_raw_Hymenoptera_limite0.05,sum0.05)->s2_raw_Hymenoptera_limite0.05
}

##**transform in present/absence table**
s2_raw_Hymenoptera_limite0.05->s2_raw_Hymenoptera_limite0.05
s2_raw_Hymenoptera_limite0.05[s2_raw_Hymenoptera_limite0.05>1]<-1 ##transform in present/absence table 

##**checking if there is any row with no presence**
s2_raw_Hymenoptera_limite0.05[,1:27]->data0.05
rowSums(data0.05)
length(which(rowSums(data0.05)!=0))
length(which(rowSums(data0.05)==0))

##**Community matrixes (samples in rows and h in cols).**
##**Hymenoptera**
t(s2_raw_Hymenoptera_limite0.05)->t_s2_f4_Hymenoptera_limite0.05 ##trasp
t_s2_f4_Hymenoptera_limite0.05[1:27,]->community_Hymenoptera_limite0.05 #NOTA_Nancy: Este numero es importante. Colocar exactamente el numero de "s2_f4[,2:52]->data0.05".
colnames(t_s2_f4_Hymenoptera_limite0.05)<-community_Hymenoptera_limite0.05[1,] 
as.data.frame(community_Hymenoptera_limite0.05)->community_Hymenoptera0.05 ##trasp including col and row names
####community_Acari[-49,]->community_Hymenoptera ##removing neg
dim(community_Hymenoptera0.05)
community_Hymenoptera0.05[order(row.names(community_Hymenoptera0.05)),]->community_Hymenoptera0.05 ##order samples
write.table (community_Hymenoptera0.05, file="../../genetic/Data_out/Hymenoptera/Hymenoptera5P/community_Hymenoptera0.05.txt") ##this is necessary for the format, not able to solve in other way
read.table ("../../genetic/Data_out/Hymenoptera/Hymenoptera5P/community_Hymenoptera0.05.txt")->community_Hymenoptera0.05


##**submatrixes by Traps in Nevado Toluca**
dim(community_Hymenoptera0.05)
community_Hymenoptera0.05[which(str_extract (row.names(community_Hymenoptera0.05), "_T_") %in% "_T_"),]->community_Hymenoptera_Traps0.05
dim(community_Hymenoptera_Traps0.05)
community_Hymenoptera_Traps0.05[,which(colSums(community_Hymenoptera_Traps0.05)!=0)]->community_Hymenoptera_Traps0.05 ##to remove no data colums
dim(community_Hymenoptera_Traps0.05)
write.table (community_Hymenoptera_Traps0.05, file="../../genetic/Data_out/Hymenoptera/Hymenoptera5P/community_Hymenoptera_Traps0.05.txt") ##this is necessary for the format, not able to solve in other way
read.table ("../../genetic/Data_out/Hymenoptera/Hymenoptera5P/community_Hymenoptera_Traps0.05.txt")->community_Hymenoptera_Traps0.05


##**Generating a general table with names and habitat parameters.**
##BY Traps
##**Generating a general table with names and habitat parameters**
row.names(community_Hymenoptera_Traps0.05)->sample_names_Mountain1_0.05
as.data.frame(sample_names_Mountain1_0.05)->sample_names_Mountain1_0.05
sample_names_Mountain1_0.05 %>% separate(sample_names_Mountain1_0.05, c("Conservation1","Mountain1","Traps","ID"), sep="_",remove=FALSE)->general_sample_Mountain1Traps0.05
general_sample_Mountain1Traps0.05
general_sample_Mountain1Traps0.05 %>% unite(Mountain1andTraps, Mountain1,Traps, sep="_",remove=FALSE)->general_sample_Mountain1Traps0.05 ####generating a variable combining layer and habitat
general_sample_Mountain1Traps0.05
write.table(general_sample_Mountain1Traps0.05, file="../../genetic/Data_out/Hymenoptera/Hymenoptera5P/general_sample_Mountain1Traps0.05.txt") ##this is the only way I found to be able to work later
read.table("../../genetic/Data_out/Hymenoptera/Hymenoptera5P/general_sample_Mountain1Traps0.05.txt",header=TRUE)->general_sample_Mountain1Traps0.05



#**HAPLOTYPE RICHNESS TABLES, PLOTS AND ANALYSES by TrapsS**
##**Hymenoptera**
as.matrix(community_Hymenoptera_Traps0.05)->community_Hymenoptera_Traps0.05
row.names(community_Hymenoptera_Traps0.05)->sample_names_Traps0.05
dim(community_Hymenoptera_Traps0.05)->dims_Traps0.05
dims_Traps0.05
.rowSums (community_Hymenoptera_Traps0.05,dims_Traps0.05[1],dims_Traps0.05[2])->sample_richness_Traps0.05 ##summatory by rows
rbind(sample_names_Traps0.05,sample_richness_Traps0.05)->richness_Traps0.05
t(richness_Traps0.05)->richness_Traps0.05
colnames(richness_Traps0.05)<-c("sample_names_Traps0.05","sample_richness_Traps0.05")
richness_Traps0.05
as.data.frame(richness_Traps0.05)->richness_Traps0.05

##**Generating variables with Traps. En la tabla rishnees. Generar variable de montana y sitio.** 
richness_Traps0.05 %>% separate(sample_names_Traps0.05, c("Conservation1","Mountain1","Traps","ID"), sep="_",remove=FALSE)->richness_Traps0.05
richness_Traps0.05
richness_Traps0.05 %>% unite(Mountain1Traps, Mountain1, Traps, sep="_",remove=FALSE)->richness_Traps0.05 ##generating a variable combining layer and habitat
richness_Traps0.05

##**Generating variables with in total Traps_C.**
#richness_Traps0.05 %>% separate(sample_names_Traps0.05, c("Conservation","Mountain1","Traps","ID"), sep="_",remove=FALSE)->richness_TrapsC
#richness_TrapsC
#richness_TrapsC %>% unite(ConservationMountain1, Conservation, Mountain1, sep="_",remove=FALSE)->richness_TrapsC ##generating a variable combining layer and habitat
#richness_TrapsC

##**BY Traps**
write.table(richness_Traps0.05, file="../../genetic/Data_out/Hymenoptera/Hymenoptera5P/richness_Traps0.05_Hymenoptera.txt") ##this is the only way I found to be able to work later
read.table("../../genetic/Data_out/Hymenoptera/Hymenoptera5P/richness_Traps0.05_Hymenoptera.txt",header=TRUE)->richness_Traps0.05

##**BY Traps_C general**
#write.table(richness_TrapsC, file="../../genetic/Data_out/Hymenoptera/Hymenoptera5P/richness_TrapsC_Hymenoptera.txt") ##this is the only way I found to be able to work later
#read.table("../../genetic/Data_out/Hymenoptera/Hymenoptera5P/richness_TrapsC_Hymenoptera.txt",header=TRUE)->richness_TrapsC

##**General plot of richness by sample in Traps**
barplot(richness_Traps0.05$sample_richness_Traps0.05,col=richness_Traps0.05$Mountain1Traps,names.arg= richness_Traps0.05$sample_names_Traps0.05,las=2,cex.names=0.5, ylab="richness_Traps0.05", main="H richness_Traps Hymenoptera_0.05")
richness_Traps0.05 %>% group_by(Mountain1Traps) %>% summarise(mean(sample_richness_Traps0.05))

##**min, max, ds Summarise**
richness_Traps0.05 %>% group_by(Mountain1Traps) %>% summarise(min(sample_richness_Traps0.05))
richness_Traps0.05 %>% group_by(Mountain1Traps) %>% summarise(max(sample_richness_Traps0.05))
richness_Traps0.05 %>% group_by(Mountain1Traps) %>% summarise(sd(sample_richness_Traps0.05))

##**General mean, min, max, ds by sample in Traps_richness_TrapsC**
#richness_TrapsC %>% group_by(ConservationMountain1) %>% summarise(mean(sample_richness_Traps0.05))
#richness_TrapsC %>% group_by(ConservationMountain1) %>% summarise(min(sample_richness_Traps0.05))
#richness_TrapsC %>% group_by(ConservationMountain1) %>% summarise(max(sample_richness_Traps0.05))
#richness_TrapsC %>% group_by(ConservationMountain1) %>% summarise(sd(sample_richness_Traps0.05))

##**Global richness by Traps.**  
plot(richness_Traps0.05$Mountain1Traps,richness_Traps0.05$sample_richness_Traps0.05,ylab="richness_Traps0.05", ylim=c(0,10), cex=1.4, cex.axis=2.3, lwd=2.5)
kruskal.test(sample_richness_Traps0.05 ~ Mountain1Traps, data = richness_Traps0.05)
#posthoc.kruskal.nemenyi.test(x=richness_Traps0.05$sample_richness_Traps0.05, g=richness_Traps0.05$Mountain1Traps, method="Bonferroni")
PMCMRplus::kwAllPairsConoverTest(x=richness_Traps0.05$sample_richness_Traps0.05, g=richness_Traps0.05$Mountain1Traps, p.adjust.method = "single-step")
##Comparison of each group against. 
#text(x=c(1,2,3,4.1), y=(49), labels=c("ab","a","ab","b"), cex=1.4)
#text(x=c(1,2,3,4), y=(2), labels=c("ab","a","ab","b"), cex=1.5)
text(x=4.5, y=0.5, labels="ns", cex=2)
mtext(c("lineages 5%"), side = 3, col = "black", line = 1, cex = 2)



#**BETADIVERSITY ORDINATIONS by Trapss** 
#**Hymenoptera**
##beta general

beta.multi(community_Hymenoptera_Traps0.05, index.family="sorensen")

##turnover by pairs, nmds, anosim
beta.pair(community_Hymenoptera_Traps0.05, index.family="sorensen")->beta.pair  ##betadiversity by pair of communities using sorensen on the precense/absence data, with estimation of turnover and nestedness datamatrixes simultaneously
metaMDS (beta.pair$beta.sim)->MDSbetasim_T0.05 ##NMDS

plot (MDSbetasim_T0.05, main="Hymenoptera_Traps_0.05") 
x<- MDSbetasim_T0.05$points[,1]
y<- MDSbetasim_T0.05$points[,2]
text(x, y, pos = 1, cex=0.7, labels = row.names (community_Hymenoptera_Traps0.05))

plot (MDSbetasim_T0.05, main="Hymenoptera_Traps_0.05")
with(general_sample_Mountain1Traps0.05,ordispider(MDSbetasim_T0.05, Traps, label=T, col="blue"))

plot (MDSbetasim_T0.05, xlim=c(-0.4, 0.4), ylim=c(-0.4, 0.4), cex.axis=1.4, cex=1.2, cex.lab=1.4, main="Hymenoptera_Traps_0.05")
with(general_sample_Mountain1Traps0.05,ordispider(MDSbetasim_T0.05, Traps, label=T, cex.lab=0.9, col= c("#153a7b", "#eaa22f", "#97518b", "#81b9d0")))

##**Anosim**
anosim(beta.pair$beta.sim, general_sample_Mountain1Traps0.05$Traps, permutations=999)

plot (MDSbetasim_T0.05, xlim=c(-0.4, 0.4), ylim=c(-0.4, 0.4), cex=2, cex.lab=1, cex.axis=2.3, lwd=4.5)
with(general_sample_Mountain1Traps0.05,ordispider(MDSbetasim_T0.05, Traps, label=T, cex.lab=1, col= c("#153a7b", "#eaa22f", "#97518b", "#81b9d0"), lwd=4.5))
##I put letter "r" in cursive and r2 value
mylabel = bquote(italic(r)^2)
text(x=0.23, y=-0.35, labels = mylabel, cex=1.5)
text(x=0.48, y=-0.35, labels="=0.14 ns", cex=1.5)
mtext(c("lineages 5%"), side = 3, col = "black", line = 1, cex = 2)


#**END TRAPS**


###################################
##**submatrixes by Groups in Nevado Toluca**
dim(community_Hymenoptera0.05)
community_Hymenoptera0.05[which(str_extract (row.names(community_Hymenoptera0.05), "G_") %in% "G_"),]->community_Hymenoptera_Groups0.05
dim(community_Hymenoptera_Groups0.05)
community_Hymenoptera_Groups0.05[,which(colSums(community_Hymenoptera_Groups0.05)!=0)]->community_Hymenoptera_Groups0.05 ##to remove no data colums
dim(community_Hymenoptera_Groups0.05)
write.table (community_Hymenoptera_Groups0.05, file="../../genetic/Data_out/Hymenoptera/Hymenoptera5P/community_Hymenoptera_Groups0.05.txt") ##this is necessary for the format, not able to solve in other way
read.table ("../../genetic/Data_out/Hymenoptera/Hymenoptera5P/community_Hymenoptera_Groups0.05.txt")->community_Hymenoptera_Groups0.05


##**Generating a general table with names and habitat parameters.**
##BY Groups
##**Generating a general table with names and habitat parameters**
row.names(community_Hymenoptera_Groups0.05)->sample_names_Mountain2_0.05
as.data.frame(sample_names_Mountain2_0.05)->sample_names_Mountain2_0.05
sample_names_Mountain2_0.05 %>% separate(sample_names_Mountain2_0.05, c("Conservation2","Mountain2","Groups","ID"), sep="_",remove=FALSE)->general_sample_Mountain2Groups0.05
general_sample_Mountain2Groups0.05
general_sample_Mountain2Groups0.05 %>% unite(Mountain2andGroups, Mountain2,Groups, sep="_",remove=FALSE)->general_sample_Mountain2Groups0.05 ####generating a variable combining layer and habitat
general_sample_Mountain2Groups0.05
write.table(general_sample_Mountain2Groups0.05, file="../../genetic/Data_out/Hymenoptera/Hymenoptera5P/general_sample_Mountain2Groups0.05.txt") ##this is the only way I found to be able to work later
read.table("../../genetic/Data_out/Hymenoptera/Hymenoptera5P/general_sample_Mountain2Groups0.05.txt",header=TRUE)->general_sample_Mountain2Groups0.05



#**HAPLOTYPE RICHNESS TABLES, PLOTS AND ANALYSES by GroupsS**
##**Hymenoptera**
as.matrix(community_Hymenoptera_Groups0.05)->community_Hymenoptera_Groups0.05
row.names(community_Hymenoptera_Groups0.05)->sample_names_Groups0.05
dim(community_Hymenoptera_Groups0.05)->dims_Groups0.05
dims_Groups0.05
.rowSums (community_Hymenoptera_Groups0.05,dims_Groups0.05[1],dims_Groups0.05[2])->sample_richness_Groups0.05 ##summatory by rows
rbind(sample_names_Groups0.05,sample_richness_Groups0.05)->richness_Groups0.05
t(richness_Groups0.05)->richness_Groups0.05
colnames(richness_Groups0.05)<-c("sample_names_Groups0.05","sample_richness_Groups0.05")
richness_Groups0.05
as.data.frame(richness_Groups0.05)->richness_Groups0.05

##**Generating variables with Groups. En la tabla rishnees. Generar variable de montana y sitio.** 
richness_Groups0.05 %>% separate(sample_names_Groups0.05, c("Conservation2","Mountain2","Groups","ID"), sep="_",remove=FALSE)->richness_Groups0.05
richness_Groups0.05
richness_Groups0.05 %>% unite(Mountain2Groups, Mountain2, Groups, sep="_",remove=FALSE)->richness_Groups0.05 ##generating a variable combining layer and habitat
richness_Groups0.05

##**Generating variables with in total Groups_C.**
#richness_Groups0.05 %>% separate(sample_names_Groups0.05, c("Conservation","Mountain1","Groups","ID"), sep="_",remove=FALSE)->richness_GroupsC
#richness_GroupsC
#richness_GroupsC %>% unite(ConservationMountain1, Conservation, Mountain1, sep="_",remove=FALSE)->richness_GroupsC ##generating a variable combining layer and habitat
#richness_GroupsC

##**BY Groups**
write.table(richness_Groups0.05, file="../../genetic/Data_out/Hymenoptera/Hymenoptera5P/richness_Groups0.05_Hymenoptera.txt") ##this is the only way I found to be able to work later
read.table("../../genetic/Data_out/Hymenoptera/Hymenoptera5P/richness_Groups0.05_Hymenoptera.txt",header=TRUE)->richness_Groups0.05

##**BY Groups_C general**
#write.table(richness_GroupsC, file="../../genetic/Data_out/Hymenoptera/Hymenoptera5P/richness_GroupsC_Hymenoptera.txt") ##this is the only way I found to be able to work later
#read.table("../../genetic/Data_out/Hymenoptera/Hymenoptera5P/richness_GroupsC_Hymenoptera.txt",header=TRUE)->richness_GroupsC

##**General plot of richness by sample in Groups**
barplot(richness_Groups0.05$sample_richness_Groups0.05,col=richness_Groups0.05$Mountain2Groups,names.arg= richness_Groups0.05$sample_names_Groups0.05,las=2,cex.names=0.5, ylab="richness_Groups0.05", main="H richness_Groups Hymenoptera_0.05")
richness_Groups0.05 %>% group_by(Mountain2Groups) %>% summarise(mean(sample_richness_Groups0.05))

##**min, max, ds Summarise**
richness_Groups0.05 %>% group_by(Mountain2Groups) %>% summarise(min(sample_richness_Groups0.05))
richness_Groups0.05 %>% group_by(Mountain2Groups) %>% summarise(max(sample_richness_Groups0.05))
richness_Groups0.05 %>% group_by(Mountain2Groups) %>% summarise(sd(sample_richness_Groups0.05))

##**General mean, min, max, ds by sample in Groups_richness_GroupsC**
#richness_GroupsC %>% group_by(ConservationMountain1) %>% summarise(mean(sample_richness_Groups0.05))
#richness_GroupsC %>% group_by(ConservationMountain1) %>% summarise(min(sample_richness_Groups0.05))
#richness_GroupsC %>% group_by(ConservationMountain1) %>% summarise(max(sample_richness_Groups0.05))
#richness_GroupsC %>% group_by(ConservationMountain1) %>% summarise(sd(sample_richness_Groups0.05))


##**Global richness by Groups.**  
plot(richness_Groups0.05$Mountain2Groups,richness_Groups0.05$sample_richness_Groups0.05,ylab="richness_Groups_0.05_Hymenoptera", cex=1.4, cex.axis=2.3, lwd=2.5, ylim=c(0,11),col= c("#153a7b", "#eaa22f", "#81b9d0"))
kruskal.test(sample_richness_Groups0.05 ~ Mountain2Groups, data = richness_Groups0.05)
#PMCMRplus::kwAllPairsNemenyiTest(x=richness_Groups_h_Hymenoptera$sample_richness_Groups_h, g=richness_Groups_h_Hymenoptera$Mountain2Groups, p.adjust.method = "single-step")
PMCMRplus::kwAllPairsConoverTest(x=richness_Groups0.05$sample_richness_Groups0.05, g=richness_Groups0.05$Mountain2Groups, p.adjust.method = "single-step")
#PMCMRplus::kwAllPairsDunnTest(x=richness_Groups_h_Hymenoptera$sample_richness_Groups_h, g=richness_Groups_h_Hymenoptera$Mountain2Groups, p.adjust.method="bonferroni")
#kwAllPairsNemenyiTest(x=richness_Groups_h_Hymenoptera$sample_richness_Groups_h, g=richness_Groups_h_Hymenoptera$Mountain2Groups)
##Comparison of each group against. 
#text(x=c(1,2,3,4,5), y=1.5, labels=c("ab","a","ab","ab","b"), cex=1.5)
text(x=3.4, y=0.5, labels="ns", cex=2)
mtext(c("lineages 5%"), side = 3, col = "black", line = 1, cex = 2)

#**BETADIVERSITY ORDINATIONS by Groupss** 
#**Hymenoptera**
##beta general

beta.multi(community_Hymenoptera_Groups0.05, index.family="sorensen")

##turnover by pairs, nmds, anosim
beta.pair(community_Hymenoptera_Groups0.05, index.family="sorensen")->beta.pair  ##betadiversity by pair of communities using sorensen on the precense/absence data, with estimation of turnover and nestedness datamatrixes simultaneously
metaMDS (beta.pair$beta.sim)->MDSbetasim_G0.05 ##NMDS

plot (MDSbetasim_G0.05, main="Hymenoptera_Groups_0.05") 
x<- MDSbetasim_G0.05$points[,1]
y<- MDSbetasim_G0.05$points[,2]
text(x, y, pos = 1, cex=0.7, labels = row.names (community_Hymenoptera_Groups0.05))

plot (MDSbetasim_G0.05, main="Hymenoptera_Groups_0.05")
with(general_sample_Mountain2Groups0.05,ordispider(MDSbetasim_G0.05, Mountain2, label=T, col="blue"))

plot (MDSbetasim_G0.05, xlim=c(-0.5, 0.5), ylim=c(-0.4, 0.4), cex.axis=1.4, cex=1.2, cex.lab=1.4, main="Hymenoptera_Groups_0.05")
with(general_sample_Mountain2Groups0.05,ordispider(MDSbetasim_G0.05, Mountain2, label=T, cex.lab=0.9, col= c("#153a7b", "#eaa22f", "#81b9d0")))

##**Anosim**
anosim(beta.pair$beta.sim, general_sample_Mountain2Groups0.05$Mountain2, permutations=999)

plot (MDSbetasim_G0.05, xlim=c(-0.4, 0.4), ylim=c(-0.4, 0.4), cex=2, cex.lab=1, cex.axis=2.3, lwd=4.5)
with(general_sample_Mountain2Groups0.05,ordispider(MDSbetasim_G0.05, Mountain2, label=T, cex.lab=1, col= c("#153a7b", "#eaa22f", "#81b9d0"), lwd=4.5))
##I put letter "r" in cursive and r2 value
mylabel = bquote(italic(r)^2)
text(x=0.23, y=-0.35, labels = mylabel, cex=1.5)
text(x=0.48, y=-0.35, labels="=0.027 ns", cex=1.5)
mtext(c("lineages 5%"), side = 3, col = "black", line = 1, cex = 2)

#**END**
