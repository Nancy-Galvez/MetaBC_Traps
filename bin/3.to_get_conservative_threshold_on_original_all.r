#**This script gets conservative thershold on origibale table OTUs/ZOTUs by artropods order.** 
###We used 6 groups: Arachnida, Coleoptera, Collembola, Diptera, Hymenoptera and Hemiptera. 
###It is important to place the correct number of columns and rows, because this can change in each group.
###In this step, we can join all the groups or put them separately as in this case.

library(dplyr)
library(tidyr)
library (knitr)

#**to_get_conservative_threshold_on_original_table_Arachnida**

##**open table with names including Region and habitat parameters**
s2_raw_all_Arachnida <- read.table("../genetic/Data_in/Arachnida/3Arachnida_138_33SOPAS_alignmentAlleleTable&SpeciesLimites.txt", sep = ",", header=TRUE)
dim(s2_raw_all_Arachnida)

##**remove additional columns and leave only names (of haplotipes), samples and taxa**
s2_raw_all_Arachnida[,c(1:31)]->s2_raw_Arachnida
dim(s2_raw_Arachnida) #'51 samples = 51 plus 1 neg (the second neg from DOM_REPS is not there because all 0)
colnames(s2_raw_Arachnida)

##**delete the ocurrences with less than 4 reads by library (same criteria than denoising)**
s2_raw_Arachnida->s2_f4_with_abundance_Arachnida 
s2_f4_with_abundance_Arachnida[s2_f4_with_abundance_Arachnida<4]<-0  #'2 warning corresponding wiht the columms of the names, taxa

##**checking if there is any row with no presence**
s2_f4_with_abundance_Arachnida[,2:31]->data_Arachnida
rowSums(data_Arachnida)
length(which(rowSums(data_Arachnida)!=0))
length(which(rowSums(data_Arachnida)==0))

##generating subtables for Acari, Arachnida, Coleoptera
###s2_f4_with_abundance[which(s2_f4_with_abundance$taxa == "Acari"),]->s2_f4_Acari 
###s2_f4_Acari[,1:50]->s2_f4_Acari #'remove taxa col
###dim(s2_f4_Acari)

###s2_f4_with_abundance[which(s2_f4_with_abundance$taxa == "Arachnida"),]->s2_f4_Arachnida 
###s2_f4_Arachnida[,1:50]->s2_f4_Arachnida #'remove taxonomy
###dim(s2_f4_Arachnida)

###s2_f4_with_abundance[which(s2_f4_with_abundance$taxa == "Coleoptera"),]->s2_f4_Coleoptera 
###s2_f4_Coleoptera[,1:50]->s2_f4_Coleoptera #'remove taxonomy
###dim(s2_f4_Coleoptera)

colSums(s2_f4_with_abundance_Arachnida[,-1])
###colSums(s2_f4_Arachnida[,-1])
###colSums(s2_f4_Coleoptera[,-1])

#**Arachnida**
##Filtrado 1% por Biblioteca
rownames(s2_f4_with_abundance_Arachnida)<-s2_f4_with_abundance_Arachnida[,1]
s2_f4_with_abundance_Arachnida[,-1]->s2_f4_Arachnida_limpia

colSums(s2_f4_Arachnida_limpia)

  for (i in 1:length(colnames(s2_f4_Arachnida_limpia)))
  {
   0.01*(sum(s2_f4_Arachnida_limpia[,i]))->lim_Arachnida
    s2_f4_Arachnida_limpia[,i][s2_f4_Arachnida_limpia[i]<lim_Arachnida]<-0
     }

colSums(s2_f4_Arachnida_limpia)

deleted_Arachnida<-which(rowSums(s2_f4_Arachnida_limpia)==0)
maintained_Arachnida<-which(rowSums(s2_f4_Arachnida_limpia)!=0)

##**Generating a table with 0 and 1 for deleted and maintained**    
names(deleted_Arachnida)->deleted_Arachnida
rep(0,length(deleted_Arachnida))->haplo_deleted_Arachnida

dim(deleted_Arachnida)<-c(length(deleted_Arachnida),1)
dim(haplo_deleted_Arachnida)<-c(length(haplo_deleted_Arachnida),1)
cbind (deleted_Arachnida, haplo_deleted_Arachnida)->haplo_deleted_Arachnida

names(maintained_Arachnida)->maintained_Arachnida  
rep(1,length(maintained_Arachnida))->haplo_maintained_Arachnida

dim(maintained_Arachnida)<-c(length(maintained_Arachnida),1)
dim(haplo_maintained_Arachnida)<-c(length(haplo_maintained_Arachnida),1)
cbind (maintained_Arachnida, haplo_maintained_Arachnida)->haplo_maintained_Arachnida

rbind(haplo_deleted_Arachnida,haplo_maintained_Arachnida)->conservative_threshold_Arachnida
table(conservative_threshold_Arachnida[,2])

rbind(conservative_threshold_Arachnida)->conservative_threshold_all_Arachnida

##**Combining conservative threshold info with intial table**
s2_raw_all_Arachnida[order(s2_raw_all_Arachnida[,1]),]->s2_raw_all_ord_Arachnida
conservative_threshold_all_Arachnida[order(conservative_threshold_all_Arachnida[,1]),]->conservative_threshold_all_ord_Arachnida
colnames(conservative_threshold_all_ord_Arachnida)<-c("haplo_names_conservative_threshold_Arachnida","conservative_threshold_Arachnida")

cbind(s2_raw_all_ord_Arachnida,conservative_threshold_all_ord_Arachnida)->s2_raw_all_Arachnida_threshold
write.table(s2_raw_all_Arachnida_threshold,file="../genetic/Data_in/Arachnida/s2_raw_all_Arachnida_threshold.txt")
#

#**to_get_conservative_threshold_on_original_table_Coleoptera**

##**open table with names including Region and habitat parameters**
s2_raw_all_Coleoptera <- read.table("../genetic/Data_in/Coleoptera/3Coleoptera_71_33SOPAS_alignment_Allele_table&SpeciesLimite.txt",sep = ",", header=TRUE)
dim(s2_raw_all_Coleoptera)

##**remove additional columns and leave only names (of haplotipes), samples and taxa**
s2_raw_all_Coleoptera[,c(1:30)]->s2_raw_Coleoptera
dim(s2_raw_Coleoptera) #'51 samples = 51 plus 1 neg (the second neg from DOM_REPS is not there because all 0)
colnames(s2_raw_Coleoptera)

##**delete the ocurrences with less than 4 reads by library (same criteria than denoising)**
s2_raw_Coleoptera->s2_f4_with_abundance_Coleoptera 
s2_f4_with_abundance_Coleoptera[s2_f4_with_abundance_Coleoptera<4]<-0  #'2 warning corresponding wiht the columms of the names, taxa

##**checking if there is any row with no presence**
s2_f4_with_abundance_Coleoptera[,2:30]->data_Coleoptera
rowSums(data_Coleoptera)
length(which(rowSums(data_Coleoptera)!=0))
length(which(rowSums(data_Coleoptera)==0))

##**generating subtables for Acari, Coleoptera, Coleoptera**
###s2_f4_with_abundance[which(s2_f4_with_abundance$taxa == "Acari"),]->s2_f4_Acari 
###s2_f4_Acari[,1:50]->s2_f4_Acari #'remove taxa col
###dim(s2_f4_Acari)

###s2_f4_with_abundance[which(s2_f4_with_abundance$taxa == "Coleoptera"),]->s2_f4_Coleoptera 
###s2_f4_Coleoptera[,1:50]->s2_f4_Coleoptera #'remove taxonomy
###dim(s2_f4_Coleoptera)

###s2_f4_with_abundance[which(s2_f4_with_abundance$taxa == "Coleoptera"),]->s2_f4_Coleoptera 
###s2_f4_Coleoptera[,1:50]->s2_f4_Coleoptera #'remove taxonomy
###dim(s2_f4_Coleoptera)

colSums(s2_f4_with_abundance_Coleoptera[,-1])
###colSums(s2_f4_Coleoptera[,-1])
###colSums(s2_f4_Coleoptera[,-1])

#**Coleoptera**
##Filtrado 1% por Biblioteca
rownames(s2_f4_with_abundance_Coleoptera)<-s2_f4_with_abundance_Coleoptera[,1]
s2_f4_with_abundance_Coleoptera[,-1]->s2_f4_Coleoptera_limpia

colSums(s2_f4_Coleoptera_limpia)

  for (i in 1:length(colnames(s2_f4_Coleoptera_limpia)))
  {
   0.01*(sum(s2_f4_Coleoptera_limpia[,i]))->lim_Coleoptera
    s2_f4_Coleoptera_limpia[,i][s2_f4_Coleoptera_limpia[i]<lim_Coleoptera]<-0
     }

colSums(s2_f4_Coleoptera_limpia)

 
deleted_Coleoptera<-which(rowSums(s2_f4_Coleoptera_limpia)==0)
maintained_Coleoptera<-which(rowSums(s2_f4_Coleoptera_limpia)!=0)

##**Generating a table with 0 and 1 for deleted and maintained**   
names(deleted_Coleoptera)->deleted_Coleoptera
rep(0,length(deleted_Coleoptera))->haplo_deleted_Coleoptera

dim(deleted_Coleoptera)<-c(length(deleted_Coleoptera),1)
dim(haplo_deleted_Coleoptera)<-c(length(haplo_deleted_Coleoptera),1)
cbind (deleted_Coleoptera, haplo_deleted_Coleoptera)->haplo_deleted_Coleoptera
  
names(maintained_Coleoptera)->maintained_Coleoptera
rep(1,length(maintained_Coleoptera))->haplo_maintained_Coleoptera

dim(maintained_Coleoptera)<-c(length(maintained_Coleoptera),1)
dim(haplo_maintained_Coleoptera)<-c(length(haplo_maintained_Coleoptera),1)
cbind (maintained_Coleoptera, haplo_maintained_Coleoptera)->haplo_maintained_Coleoptera

rbind(haplo_deleted_Coleoptera,haplo_maintained_Coleoptera)->conservative_threshold_Coleoptera
table(conservative_threshold_Coleoptera[,2])

rbind(conservative_threshold_Coleoptera)->conservative_threshold_all_Coleoptera

##**Combining conservative threshold info with intial table**
s2_raw_all_Coleoptera[order(s2_raw_all_Coleoptera[,1]),]->s2_raw_all_ord_Coleoptera
conservative_threshold_all_Coleoptera[order(conservative_threshold_all_Coleoptera[,1]),]->conservative_threshold_all_ord_Coleoptera
colnames(conservative_threshold_all_ord_Coleoptera)<-c("haplo_names_conservative_threshold_Coleoptera","conservative_threshold_Coleoptera")

cbind(s2_raw_all_ord_Coleoptera,conservative_threshold_all_ord_Coleoptera)->s2_raw_all_Coleoptera_threshold
write.table(s2_raw_all_Coleoptera_threshold,file="../genetic/Data_in/Coleoptera/s2_raw_all_Coleoptera_threshold.txt")
#

#**to_get_conservative_threshold_on_original_table_Collembola**

##**open table with names including Region and habitat parameters**
s2_raw_all_Collembola <- read.table("../genetic/Data_in/Collembola/3Collembola_158_33SOPAS_alignment_allele_table&SpeciesLimite.txt",sep = ",", header=TRUE)
dim(s2_raw_all_Collembola)

##**remove additional columns and leave only names (of haplotipes), samples and taxa**
s2_raw_all_Collembola[,c(1:26)]->s2_raw_Collembola
dim(s2_raw_Collembola) #'51 samples = 51 plus 1 neg (the second neg from DOM_REPS is not there because all 0)
colnames(s2_raw_Collembola)

##**delete the ocurrences with less than 4 reads by library (same criteria than denoising)**
s2_raw_Collembola->s2_f4_with_abundance_Collembola 
s2_f4_with_abundance_Collembola[s2_f4_with_abundance_Collembola<4]<-0  #'2 warning corresponding wiht the columms of the names, taxa

##**checking if there is any row with no presence**
s2_f4_with_abundance_Collembola[,2:26]->data_Collembola
rowSums(data_Collembola)
length(which(rowSums(data_Collembola)!=0))
length(which(rowSums(data_Collembola)==0))

##**generating subtables for Acari, Collembola, Coleoptera**
###s2_f4_with_abundance[which(s2_f4_with_abundance$taxa == "Acari"),]->s2_f4_Acari 
###s2_f4_Acari[,1:50]->s2_f4_Acari #'remove taxa col
###dim(s2_f4_Acari)

###s2_f4_with_abundance[which(s2_f4_with_abundance$taxa == "Collembola"),]->s2_f4_Collembola 
###s2_f4_Collembola[,1:50]->s2_f4_Collembola #'remove taxonomy
###dim(s2_f4_Collembola)

###s2_f4_with_abundance[which(s2_f4_with_abundance$taxa == "Coleoptera"),]->s2_f4_Coleoptera 
###s2_f4_Coleoptera[,1:50]->s2_f4_Coleoptera #'remove taxonomy
###dim(s2_f4_Coleoptera)

colSums(s2_f4_with_abundance_Collembola[,-1])
###colSums(s2_f4_Collembola[,-1])
###colSums(s2_f4_Coleoptera[,-1])

#**Collembola**
###Filtrado 1% por Biblioteca
rownames(s2_f4_with_abundance_Collembola)<-s2_f4_with_abundance_Collembola[,1]
s2_f4_with_abundance_Collembola[,-1]->s2_f4_Collembola_limpia

colSums(s2_f4_Collembola_limpia)

  for (i in 1:length(colnames(s2_f4_Collembola_limpia)))
  {
   0.01*(sum(s2_f4_Collembola_limpia[,i]))->lim_Collembola
    s2_f4_Collembola_limpia[,i][s2_f4_Collembola_limpia[i]<lim_Collembola]<-0
     }

colSums(s2_f4_Collembola_limpia)
 
deleted_Collembola<-which(rowSums(s2_f4_Collembola_limpia)==0)
maintained_Collembola<-which(rowSums(s2_f4_Collembola_limpia)!=0)

##**Generating a table with 0 and 1 for deleted and maintained**    
names(deleted_Collembola)->deleted_Collembola
rep(0,length(deleted_Collembola))->haplo_deleted_Collembola

dim(deleted_Collembola)<-c(length(deleted_Collembola),1)
dim(haplo_deleted_Collembola)<-c(length(haplo_deleted_Collembola),1)
cbind (deleted_Collembola, haplo_deleted_Collembola)->haplo_deleted_Collembola
  
names(maintained_Collembola)->maintained_Collembola
rep(1,length(maintained_Collembola))->haplo_maintained_Collembola

dim(maintained_Collembola)<-c(length(maintained_Collembola),1)
dim(haplo_maintained_Collembola)<-c(length(haplo_maintained_Collembola),1)
cbind (maintained_Collembola, haplo_maintained_Collembola)->haplo_maintained_Collembola

rbind(haplo_deleted_Collembola,haplo_maintained_Collembola)->conservative_threshold_Collembola
table(conservative_threshold_Collembola[,2])

rbind(conservative_threshold_Collembola)->conservative_threshold_all_Collembola

##**Combining conservative threshold info with intial table**
s2_raw_all_Collembola[order(s2_raw_all_Collembola[,1]),]->s2_raw_all_ord_Collembola
conservative_threshold_all_Collembola[order(conservative_threshold_all_Collembola[,1]),]->conservative_threshold_all_ord_Collembola
colnames(conservative_threshold_all_ord_Collembola)<-c("haplo_names_conservative_threshold_Collembola","conservative_threshold_Collembola")

cbind(s2_raw_all_ord_Collembola,conservative_threshold_all_ord_Collembola)->s2_raw_all_Collembola_threshold
write.table(s2_raw_all_Collembola_threshold,file="../genetic/Data_in/Collembola/s2_raw_all_Collembola_threshold.txt")
#

#**to_get_conservative_threshold_on_original_table_Diptera**

##**open table with names including Region and habitat parameters**
s2_raw_all_Diptera <- read.table("../genetic/Data_in/Diptera/3Dipteras_583_33SOPAS_alignment_Allele_table&SpeciesLimite.txt", sep = ",", header=TRUE)
dim(s2_raw_all_Diptera)

##**remove additional columns and leave only names (of haplotipes), samples and taxa**
s2_raw_all_Diptera[,c(1:34)]->s2_raw_Diptera
dim(s2_raw_Diptera) #'51 samples = 51 plus 1 neg (the second neg from DOM_REPS is not there because all 0)
colnames(s2_raw_Diptera)

##**delete the ocurrences with less than 4 reads by library (same criteria than denoising)**
s2_raw_Diptera->s2_f4_with_abundance_Diptera
s2_f4_with_abundance_Diptera[s2_f4_with_abundance_Diptera<4]<-0  #'2 warning corresponding wiht the columms of the names, taxa

##**checking if there is any row with no presence**
s2_f4_with_abundance_Diptera[,2:34]->data_Diptera
rowSums(data_Diptera)
length(which(rowSums(data_Diptera)!=0))
length(which(rowSums(data_Diptera)==0))

##**generating subtables for Acari, Diptera, Coleoptera**
###s2_f4_with_abundance[which(s2_f4_with_abundance$taxa == "Acari"),]->s2_f4_Acari 
###s2_f4_Acari[,1:50]->s2_f4_Acari #'remove taxa col
###dim(s2_f4_Acari)

###s2_f4_with_abundance[which(s2_f4_with_abundance$taxa == "Diptera"),]->s2_f4_Diptera 
###s2_f4_Diptera[,1:50]->s2_f4_Diptera #'remove taxonomy
###dim(s2_f4_Diptera)

###s2_f4_with_abundance[which(s2_f4_with_abundance$taxa == "Coleoptera"),]->s2_f4_Coleoptera 
###s2_f4_Coleoptera[,1:50]->s2_f4_Coleoptera #'remove taxonomy
###dim(s2_f4_Coleoptera)

colSums(s2_f4_with_abundance_Diptera[,-1])
#colSums(s2_f4_Diptera[,-1])
#colSums(s2_f4_Coleoptera[,-1])

#Diptera
#Filtrado 1% por Biblioteca

rownames(s2_f4_with_abundance_Diptera)<-s2_f4_with_abundance_Diptera[,1]
s2_f4_with_abundance_Diptera[,-1]->s2_f4_Diptera_limpia

colSums(s2_f4_Diptera_limpia)

  for (i in 1:length(colnames(s2_f4_Diptera_limpia)))
  {
   0.01*(sum(s2_f4_Diptera_limpia[,i]))->lim_Diptera
    s2_f4_Diptera_limpia[,i][s2_f4_Diptera_limpia[i]<lim_Diptera]<-0
     }

colSums(s2_f4_Diptera_limpia)
 
deleted_Diptera<-which(rowSums(s2_f4_Diptera_limpia)==0)
maintained_Diptera<-which(rowSums(s2_f4_Diptera_limpia)!=0)

#Generating a table with 0 and 1 for deleted and maintained    
names(deleted_Diptera)->deleted_Diptera
rep(0,length(deleted_Diptera))->haplo_deleted_Diptera

dim(deleted_Diptera)<-c(length(deleted_Diptera),1)
dim(haplo_deleted_Diptera)<-c(length(haplo_deleted_Diptera),1)
cbind (deleted_Diptera, haplo_deleted_Diptera)->haplo_deleted_Diptera

names(maintained_Diptera)->maintained_Diptera
rep(1,length(maintained_Diptera))->haplo_maintained_Diptera

dim(maintained_Diptera)<-c(length(maintained_Diptera),1)
dim(haplo_maintained_Diptera)<-c(length(haplo_maintained_Diptera),1)
cbind (maintained_Diptera, haplo_maintained_Diptera)->haplo_maintained_Diptera

rbind(haplo_deleted_Diptera,haplo_maintained_Diptera)->conservative_threshold_Diptera
table(conservative_threshold_Diptera[,2])

rbind(conservative_threshold_Diptera)->conservative_threshold_all_Diptera

##**Combining conservative threshold info with intial table**
s2_raw_all_Diptera[order(s2_raw_all_Diptera[,1]),]->s2_raw_all_ord_Diptera
conservative_threshold_all_Diptera[order(conservative_threshold_all_Diptera[,1]),]->conservative_threshold_all_ord_Diptera
colnames(conservative_threshold_all_ord_Diptera)<-c("haplo_names_conservative_threshold_Diptera","conservative_threshold_Diptera")

cbind(s2_raw_all_ord_Diptera,conservative_threshold_all_ord_Diptera)->s2_raw_all_Diptera_threshold
write.table(s2_raw_all_Diptera_threshold,file="../genetic/Data_in/Diptera/s2_raw_all_Diptera_threshold.txt")
#

#**to_get_conservative_threshold_on_original_table_Hemiptera**

##**open table with names including Region and habitat parameters**
s2_raw_all_Hemiptera <- read.table("../genetic/Data_in/Hemiptera/3Hemiptera_150_33SOPAS_alignment_allele_table&SpeciesLimite.txt", sep = ",", header=TRUE)
dim(s2_raw_all_Hemiptera)

##**remove additional columns and leave only names (of haplotipes), samples and taxa**
s2_raw_all_Hemiptera[,c(1:29)]->s2_raw_Hemiptera
dim(s2_raw_Hemiptera) #'51 samples = 51 plus 1 neg (the second neg from DOM_REPS is not there because all 0)
colnames(s2_raw_Hemiptera)

##**delete the ocurrences with less than 4 reads by library (same criteria than denoising)**
s2_raw_Hemiptera->s2_f4_with_abundance_Hemiptera
s2_f4_with_abundance_Hemiptera[s2_f4_with_abundance_Hemiptera<4]<-0  #'2 warning corresponding wiht the columms of the names, taxa

##**checking if there is any row with no presence**
s2_f4_with_abundance_Hemiptera[,2:29]->data_Hemiptera
rowSums(data_Hemiptera)
length(which(rowSums(data_Hemiptera)!=0))
length(which(rowSums(data_Hemiptera)==0))

##**generating subtables for Acari, Hemiptera, Coleoptera**
###s2_f4_with_abundance[which(s2_f4_with_abundance$taxa == "Acari"),]->s2_f4_Acari 
###s2_f4_Acari[,1:50]->s2_f4_Acari #'remove taxa col
###dim(s2_f4_Acari)

###s2_f4_with_abundance[which(s2_f4_with_abundance$taxa == "Hemiptera"),]->s2_f4_Hemiptera 
###s2_f4_Hemiptera[,1:50]->s2_f4_Hemiptera #'remove taxonomy
###dim(s2_f4_Hemiptera)

###s2_f4_with_abundance[which(s2_f4_with_abundance$taxa == "Coleoptera"),]->s2_f4_Coleoptera 
###s2_f4_Coleoptera[,1:50]->s2_f4_Coleoptera #'remove taxonomy
###dim(s2_f4_Coleoptera)

colSums(s2_f4_with_abundance_Hemiptera[,-1])
#colSums(s2_f4_Hemiptera[,-1])
#colSums(s2_f4_Coleoptera[,-1])

#**Hemiptera**
##Filtrado 1% por Biblioteca
rownames(s2_f4_with_abundance_Hemiptera)<-s2_f4_with_abundance_Hemiptera[,1]
s2_f4_with_abundance_Hemiptera[,-1]->s2_f4_Hemiptera_limpia

colSums(s2_f4_Hemiptera_limpia)

  for (i in 1:length(colnames(s2_f4_Hemiptera_limpia)))
  {
   0.01*(sum(s2_f4_Hemiptera_limpia[,i]))->lim_Hemiptera
    s2_f4_Hemiptera_limpia[,i][s2_f4_Hemiptera_limpia[i]<lim_Hemiptera]<-0
     }

colSums(s2_f4_Hemiptera_limpia)

deleted_Hemiptera<-which(rowSums(s2_f4_Hemiptera_limpia)==0)
maintained_Hemiptera<-which(rowSums(s2_f4_Hemiptera_limpia)!=0)

##**Generating a table with 0 and 1 for deleted and maintained**
names(deleted_Hemiptera)->deleted_Hemiptera
rep(0,length(deleted_Hemiptera))->haplo_deleted_Hemiptera

dim(deleted_Hemiptera)<-c(length(deleted_Hemiptera),1)
dim(haplo_deleted_Hemiptera)<-c(length(haplo_deleted_Hemiptera),1)
cbind (deleted_Hemiptera, haplo_deleted_Hemiptera)->haplo_deleted_Hemiptera

names(maintained_Hemiptera)->maintained_Hemiptera  
rep(1,length(maintained_Hemiptera))->haplo_maintained_Hemiptera

dim(maintained_Hemiptera)<-c(length(maintained_Hemiptera),1)
dim(haplo_maintained_Hemiptera)<-c(length(haplo_maintained_Hemiptera),1)
cbind (maintained_Hemiptera, haplo_maintained_Hemiptera)->haplo_maintained_Hemiptera

rbind(haplo_deleted_Hemiptera,haplo_maintained_Hemiptera)->conservative_threshold_Hemiptera
table(conservative_threshold_Hemiptera[,2])

rbind(conservative_threshold_Hemiptera)->conservative_threshold_all_Hemiptera

##**Combining conservative threshold info with intial table**
s2_raw_all_Hemiptera[order(s2_raw_all_Hemiptera[,1]),]->s2_raw_all_ord_Hemiptera
conservative_threshold_all_Hemiptera[order(conservative_threshold_all_Hemiptera[,1]),]->conservative_threshold_all_ord_Hemiptera
colnames(conservative_threshold_all_ord_Hemiptera)<-c("haplo_names_conservative_threshold_Hemiptera","conservative_threshold_Hemiptera")

cbind(s2_raw_all_ord_Hemiptera,conservative_threshold_all_ord_Hemiptera)->s2_raw_all_Hemiptera_threshold
write.table(s2_raw_all_Hemiptera_threshold,file="../genetic/Data_in/Hemiptera/s2_raw_all_Hemiptera_threshold.txt")
#

#**to_get_conservative_threshold_on_original_table_Hymenoptera**

##**open table with names including Region and habitat parameters**
s2_raw_all_Hymenoptera <- read.table("../genetic/Data_in/Hymenoptera/3Hymenoptera_161_33SOPAS_alignment_allele_table&SpeciesLimite.txt",sep = ",", header=TRUE)
dim(s2_raw_all_Hymenoptera)

##**remove additional columns and leave only names (of haplotipes), samples and taxa**
s2_raw_all_Hymenoptera[,c(1:28)]->s2_raw_Hymenoptera
dim(s2_raw_Hymenoptera) #'51 samples = 51 plus 1 neg (the second neg from DOM_REPS is not there because all 0)
colnames(s2_raw_Hymenoptera)

##**delete the ocurrences with less than 4 reads by library (same criteria than denoising)**
s2_raw_Hymenoptera->s2_f4_with_abundance_Hymenoptera 
s2_f4_with_abundance_Hymenoptera[s2_f4_with_abundance_Hymenoptera<4]<-0  #'2 warning corresponding wiht the columms of the names, taxa

##**checking if there is any row with no presence**
s2_f4_with_abundance_Hymenoptera[,2:28]->data_Hymenoptera
rowSums(data_Hymenoptera)
length(which(rowSums(data_Hymenoptera)!=0))
length(which(rowSums(data_Hymenoptera)==0))

##**generating subtables for Acari, Hymenoptera, Coleoptera**
###s2_f4_with_abundance[which(s2_f4_with_abundance$taxa == "Acari"),]->s2_f4_Acari 
###s2_f4_Acari[,1:50]->s2_f4_Acari #'remove taxa col
###dim(s2_f4_Acari)

###s2_f4_with_abundance[which(s2_f4_with_abundance$taxa == "Hymenoptera"),]->s2_f4_Hymenoptera 
###s2_f4_Hymenoptera[,1:50]->s2_f4_Hymenoptera #'remove taxonomy
###dim(s2_f4_Hymenoptera)

###s2_f4_with_abundance[which(s2_f4_with_abundance$taxa == "Coleoptera"),]->s2_f4_Coleoptera 
###s2_f4_Coleoptera[,1:50]->s2_f4_Coleoptera #'remove taxonomy
###dim(s2_f4_Coleoptera)

colSums(s2_f4_with_abundance_Hymenoptera[,-1])
###colSums(s2_f4_Hymenoptera[,-1])
###colSums(s2_f4_Coleoptera[,-1])

#**Hymenoptera**
##Filtrado 1% por Biblioteca
rownames(s2_f4_with_abundance_Hymenoptera)<-s2_f4_with_abundance_Hymenoptera[,1]
s2_f4_with_abundance_Hymenoptera[,-1]->s2_f4_Hymenoptera_limpia

colSums(s2_f4_Hymenoptera_limpia)

  for (i in 1:length(colnames(s2_f4_Hymenoptera_limpia)))
  {
   0.01*(sum(s2_f4_Hymenoptera_limpia[,i]))->lim_Hymenoptera
    s2_f4_Hymenoptera_limpia[,i][s2_f4_Hymenoptera_limpia[i]<lim_Hymenoptera]<-0
     }

colSums(s2_f4_Hymenoptera_limpia)

deleted_Hymenoptera<-which(rowSums(s2_f4_Hymenoptera_limpia)==0)
maintained_Hymenoptera<-which(rowSums(s2_f4_Hymenoptera_limpia)!=0)

##**Generating a table with 0 and 1 for deleted and maintained**    
names(deleted_Hymenoptera)->deleted_Hymenoptera
rep(0,length(deleted_Hymenoptera))->haplo_deleted_Hymenoptera

dim(deleted_Hymenoptera)<-c(length(deleted_Hymenoptera),1)
dim(haplo_deleted_Hymenoptera)<-c(length(haplo_deleted_Hymenoptera),1)
cbind (deleted_Hymenoptera, haplo_deleted_Hymenoptera)->haplo_deleted_Hymenoptera

names(maintained_Hymenoptera)->maintained_Hymenoptera
rep(1,length(maintained_Hymenoptera))->haplo_maintained_Hymenoptera

dim(maintained_Hymenoptera)<-c(length(maintained_Hymenoptera),1)
dim(haplo_maintained_Hymenoptera)<-c(length(haplo_maintained_Hymenoptera),1)
cbind (maintained_Hymenoptera, haplo_maintained_Hymenoptera)->haplo_maintained_Hymenoptera

rbind(haplo_deleted_Hymenoptera,haplo_maintained_Hymenoptera)->conservative_threshold_Hymenoptera
table(conservative_threshold_Hymenoptera[,2])

rbind(conservative_threshold_Hymenoptera)->conservative_threshold_all_Hymenoptera

##**Combining conservative threshold info with intial table**
s2_raw_all_Hymenoptera[order(s2_raw_all_Hymenoptera[,1]),]->s2_raw_all_ord_Hymenoptera
conservative_threshold_all_Hymenoptera[order(conservative_threshold_all_Hymenoptera[,1]),]->conservative_threshold_all_ord_Hymenoptera
colnames(conservative_threshold_all_ord_Hymenoptera)<-c("haplo_names_conservative_threshold_Hymenoptera","conservative_threshold_Hymenoptera")

cbind(s2_raw_all_ord_Hymenoptera,conservative_threshold_all_ord_Hymenoptera)->s2_raw_all_Hymenoptera_threshold
write.table(s2_raw_all_Hymenoptera_threshold,file="../genetic/Data_in/Hymenoptera/s2_raw_all_Hymenoptera_threshold.txt")
#

#**to_get_conservative_threshold_on_original_table_Myriapoda**

##**open table with names including Region and habitat parameters**
#s2_raw_all_Myriapoda <- read.table("../genetic/Data_in/Myriapoda/AlleleTableyLimiteEspeciesMyriapoda.txt",sep = ",", header=TRUE)
#dim(s2_raw_all_Myriapoda)

##**remove additional columns and leave only names (of haplotipes), samples and taxa**
#s2_raw_all_Myriapoda[,c(1:12)]->s2_raw_Myriapoda
#dim(s2_raw_Myriapoda) #'51 samples = 51 plus 1 neg (the second neg from DOM_REPS is not there because all 0)
#colnames(s2_raw_Myriapoda)

##**delete the ocurrences with less than 4 reads by library (same criteria than denoising)**
#s2_raw_Myriapoda->s2_f4_with_abundance_Myriapoda 
#s2_f4_with_abundance_Myriapoda[s2_f4_with_abundance_Myriapoda<4]<-0  #'2 warning corresponding wiht the columms of the names, taxa

##**checking if there is any row with no presence**
#s2_f4_with_abundance_Myriapoda[,2:12]->data_Myriapoda
#rowSums(data_Myriapoda)
#length(which(rowSums(data_Myriapoda)!=0))
#length(which(rowSums(data_Myriapoda)==0))

##**generating subtables for Acari, Myriapoda, Coleoptera**
###s2_f4_with_abundance[which(s2_f4_with_abundance$taxa == "Acari"),]->s2_f4_Acari 
###s2_f4_Acari[,1:50]->s2_f4_Acari #'remove taxa col
###dim(s2_f4_Acari)

###s2_f4_with_abundance[which(s2_f4_with_abundance$taxa == "Myriapoda"),]->s2_f4_Myriapoda 
###s2_f4_Myriapoda[,1:50]->s2_f4_Myriapoda #'remove taxonomy
###dim(s2_f4_Myriapoda)

###s2_f4_with_abundance[which(s2_f4_with_abundance$taxa == "Coleoptera"),]->s2_f4_Coleoptera 
###s2_f4_Coleoptera[,1:50]->s2_f4_Coleoptera #'remove taxonomy
###dim(s2_f4_Coleoptera)

#colSums(s2_f4_with_abundance_Myriapoda[,-1])
###colSums(s2_f4_Myriapoda[,-1])
###colSums(s2_f4_Coleoptera[,-1])

#**Myriapoda**
##Filtrado 1% por Biblioteca
#rownames(s2_f4_with_abundance_Myriapoda)<-s2_f4_with_abundance_Myriapoda[,1]
#s2_f4_with_abundance_Myriapoda[,-1]->s2_f4_Myriapoda_limpia

#colSums(s2_f4_Myriapoda_limpia)

#  for (i in 1:length(colnames(s2_f4_Myriapoda_limpia)))
#  {
#   0.01*(sum(s2_f4_Myriapoda_limpia[,i]))->lim_Myriapoda
#    s2_f4_Myriapoda_limpia[,i][s2_f4_Myriapoda_limpia[i]<lim_Myriapoda]<-0
#     }

#colSums(s2_f4_Myriapoda_limpia)
 
#deleted_Myriapoda<-which(rowSums(s2_f4_Myriapoda_limpia)==0)
#maintained_Myriapoda<-which(rowSums(s2_f4_Myriapoda_limpia)!=0)

##**Generating a table with 0 and 1 for deleted and maintained**    
#names(deleted_Myriapoda)->deleted_Myriapoda
#rep(0,length(deleted_Myriapoda))->haplo_deleted_Myriapoda

#dim(deleted_Myriapoda)<-c(length(deleted_Myriapoda),1)
#dim(haplo_deleted_Myriapoda)<-c(length(haplo_deleted_Myriapoda),1)
#cbind (deleted_Myriapoda, haplo_deleted_Myriapoda)->haplo_deleted_Myriapoda
  
#names(maintained_Myriapoda)->maintained_Myriapoda  
#rep(1,length(maintained_Myriapoda))->haplo_maintained_Myriapoda

#dim(maintained_Myriapoda)<-c(length(maintained_Myriapoda),1)
#dim(haplo_maintained_Myriapoda)<-c(length(haplo_maintained_Myriapoda),1)
#cbind (maintained_Myriapoda, haplo_maintained_Myriapoda)->haplo_maintained_Myriapoda

#rbind(haplo_deleted_Myriapoda,haplo_maintained_Myriapoda)->conservative_threshold_Myriapoda
#table(conservative_threshold_Myriapoda[,2])

#rbind(conservative_threshold_Myriapoda)->conservative_threshold_all_Myriapoda

##**Combining conservative threshold info with intial table**
#s2_raw_all_Myriapoda[order(s2_raw_all_Myriapoda[,1]),]->s2_raw_all_ord_Myriapoda
#conservative_threshold_all_Myriapoda[order(conservative_threshold_all_Myriapoda[,1]),]->conservative_threshold_all_ord_Myriapoda
#colnames(conservative_threshold_all_ord_Myriapoda)<-c("haplo_names_conservative_threshold_Myriapoda","conservative_threshold_Myriapoda")

#cbind(s2_raw_all_ord_Myriapoda,conservative_threshold_all_ord_Myriapoda)->s2_raw_all_Myriapoda_threshold
#write.table(s2_raw_all_Myriapoda_threshold,file="../genetic/Data_in/Myriapoda/s2_raw_all_Myriapoda_threshold.txt")

#**END**
