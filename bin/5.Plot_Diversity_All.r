library(raster)
require(gridExtra)
library(sp)
library(ggplot2)
library(patchwork)
library(devtools)
library(easypackages)
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

# gert rcmdcheck usethis devtools easypackages
#**Community diversity and composition at the haplotype, 3% and 5% lineages for each of the eight taxonomic orders studied**

#####This script gets all scripts of diversity using 8 arthropods order at multi-hierarchical levels 
#####We used 8 groups: Arachnda, Coleoptera, Collembola, Diptera, Hemiptera, Hymenoptera, Lepidoptera, and Myriapoda. 

#**plots Community diversity and composition**

png(filename="../figures/Figure1_Richness_Traps.png", width=820 , height=1000, units="px") # set size of the file to plot 
par(mfrow=c(5,3), mar = c(3, 2.5, 2, 2), omi=c(0.5, 0.5, 0.5, 0.5)) #the number of rows and columns the figure would have

##**Creating plot**
##**Community diversity and composition at the haplotype, 3% and 5% lineages for each of the eight taxonomic orders studied**

#####This script gets all scripts of diversity using 6 arthropods order at multi-hierarchical levels 
#####We used 6 groups: Arachnda, Coleoptera, Collembola, Diptera, Hemiptera, and Hymenoptera.

#**plot Global Richness by Sites**
#**Diptera**

read.table("../genetic/Data_out/Diptera/Diptera_Haplotypes/richness_Traps_h_Diptera.txt",header=TRUE)->richness_Traps_h_Diptera

plot(richness_Traps_h_Diptera$Mountain1Traps,richness_Traps_h_Diptera$sample_richness_Traps_h, xaxt="n", ylim=c(0,60),cex=1.4, cex.axis=2.3, lwd=2.5)
kruskal.test(sample_richness_Traps_h ~ Mountain1Traps, data = richness_Traps_h_Diptera)
PMCMRplus::kwAllPairsConoverTest(x=richness_Traps_h_Diptera$sample_richness_Traps_h, g=richness_Traps_h_Diptera$Mountain1Traps, p.adjust.method = "single-step")
##Comparison of each group against. 
axis(1, at=NULL, label=F, tick=TRUE, lwd=1.3)
text(x=c(1,2,3,4), y=(2), labels=c("ab","a","ab","b"), cex=2)
text(x=4.3, y=1, labels="*", cex=3.5)
mtext(c("haplotypes"), side = 3, col = "black", line = 1, cex=2)

read.table("../genetic/Data_out/Diptera/Diptera3P/richness_Traps0.03_Diptera.txt",header=TRUE)->richness_Traps0.03

plot(richness_Traps0.03$Mountain1Traps,richness_Traps0.03$sample_richness_Traps0.03, xaxt="n", ylim=c(0,30), cex=1.4, cex.axis=2.3, lwd=2.5)
kruskal.test(sample_richness_Traps0.03 ~ Mountain1Traps, data = richness_Traps0.03)
PMCMRplus::kwAllPairsConoverTest(x=richness_Traps0.03$sample_richness_Traps0.03, g=richness_Traps0.03$Mountain1Traps, p.adjust.method = "single-step")
##Comparison of each group against.    
axis(1, at=NULL, label=F, tick=TRUE, lwd=1.3)
text(x=c(1,2,3,4), y=(2), labels=c("ab","a","ab","b"), cex=2)
text(x=4.3, y=1, labels="*", cex=3.5)
mtext(c("3% lineages"), side = 3, col = "black", line = 1, cex=2)


read.table("../genetic/Data_out/Diptera/Diptera5P/richness_Traps0.05_Diptera.txt",header=TRUE)->richness_Traps0.05

plot(richness_Traps0.05$Mountain1Traps,richness_Traps0.05$sample_richness_Traps0.05, xaxt="n", ylim=c(0,30), cex=1.4, cex.axis=2.3, lwd=2.5)
kruskal.test(sample_richness_Traps0.05 ~ Mountain1Traps, data = richness_Traps0.05)
PMCMRplus::kwAllPairsConoverTest(x=richness_Traps0.05$sample_richness_Traps0.05, g=richness_Traps0.05$Mountain1Traps, p.adjust.method = "single-step")
##Comparison of each group against. 
#text(x=c(1,2,3,4.1), y=(49), labels=c("ab","a","ab","b"), cex=1.4)
axis(1, at=NULL, label=F, tick=TRUE, lwd=1.3)
text(x=c(1,2,3,4), y=(2), labels=c("a","a","a","b"), cex=2)
text(x=4.3, y=1, labels="*", cex=3.5)
mtext(c("5% lineages"), side = 3, col = "black", line = 1, cex=2)
mtext(c("Diptera"), side = 4, col = "black", line = 2, cex=2)

#
#**Collembola**
#read.table("../genetic/Data_out/Collembola/Collembola_Haplotypes/richness_Traps_h_Collembola.txt",header=TRUE)->richness_Traps_h_Collembola
#plot(richness_Traps_h_Collembola$Mountain1Traps,richness_Traps_h_Collembola$sample_richness_Traps_h, xaxt="n", cex=1.4, cex.axis=2.3, lwd=2.5, ylim=c(0,30))
#kruskal.test(sample_richness_Traps_h ~ Mountain1Traps, data = richness_Traps_h_Collembola)
#PMCMRplus::kwAllPairsConoverTest(x=richness_Traps_h_Collembola$sample_richness_Traps_h, g=richness_Traps_h_Collembola$Mountain1Traps, p.adjust.method = "single-step")
##Comparison of each group against. 
#text(x=c(1,2,3,4), y=(2), labels=c("ab","a","ab","b"), cex=1.5)
#axis(1, at=NULL, label=F, tick=TRUE, lwd=1.3)
#text(x=4.3, y=0.5, labels="ns", cex=2)

#read.table("../genetic/Data_out/Collembola/Collembola3P/richness_Traps0.03_Collembola.txt",header=TRUE)->richness_Traps0.03

#plot(richness_Traps0.03$Mountain1Traps,richness_Traps0.03$sample_richness_Traps0.03, xaxt="n", ylim=c(0,15), cex=1.4, cex.axis=2.3, lwd=2.5)
#kruskal.test(sample_richness_Traps0.03 ~ Mountain1Traps, data = richness_Traps0.03)
#posthoc.kruskal.nemenyi.test(x=richness_Traps0.03$sample_richness_Traps0.03, g=richness_Traps0.03$Mountain1Traps, method="Bonferroni")
#PMCMRplus::kwAllPairsConoverTest(x=richness_Traps0.03$sample_richness_Traps0.03, g=richness_Traps0.03$Mountain1Traps, p.adjust.method = "single-step")
##Comparison of each group against. 
#axis(1, at=NULL, label=F, tick=TRUE, lwd=1.3)
#text(x=c(1,2,3,4.1), y=(49), labels=c("ab","a","ab","b"), cex=1.4)
#text(x=c(1,2,3,4), y=(2), labels=c("ab","a","ab","b"), cex=1.5)
#text(x=4.3, y=0.5, labels="ns", cex=3)

#read.table("../genetic/Data_out/Collembola/Collembola5P/richness_Traps0.05_Collembola.txt",header=TRUE)->richness_Traps0.05


#plot(richness_Traps0.05$Mountain1Traps,richness_Traps0.05$sample_richness_Traps0.05, xaxt="n", ylim=c(0,15), cex=1.4, cex.axis=2.3, lwd=2.5)
#kruskal.test(sample_richness_Traps0.05 ~ Mountain1Traps, data = richness_Traps0.05)
#posthoc.kruskal.nemenyi.test(x=richness_Traps0.05$sample_richness_Traps0.05, g=richness_Traps0.05$Mountain1Traps, method="Bonferroni")
#PMCMRplus::kwAllPairsConoverTest(x=richness_Traps0.05$sample_richness_Traps0.05, g=richness_Traps0.05$Mountain1Traps, p.adjust.method = "single-step")
##Comparison of each group against. 
#axis(1, at=NULL, label=F, tick=TRUE, lwd=1.3)
#text(x=c(1,2,3,4.1), y=(49), labels=c("ab","a","ab","b"), cex=1.4)
#text(x=c(1,2,3,4), y=(2), labels=c("ab","a","ab","b"), cex=1.5)
#text(x=4.3, y=0.5, labels="ns", cex=2)
#mtext(c("Collembola"), side = 4, col = "black", line = 2, cex = 2)
#

#**Arachnida**
read.table("../genetic/Data_out/Arachnida/Arachnida_Haplotypes/richness_Traps_h_Arachnida.txt",header=TRUE)->richness_Traps_h_Arachnida


plot(richness_Traps_h_Arachnida$Mountain1Traps,richness_Traps_h_Arachnida$sample_richness_Traps_h, xaxt="n", cex=1.4, cex.axis=2.3, lwd=2.5, ylim=c(0,20))
kruskal.test(sample_richness_Traps_h ~ Mountain1Traps, data = richness_Traps_h_Arachnida)
PMCMRplus::kwAllPairsConoverTest(x=richness_Traps_h_Arachnida$sample_richness_Traps_h, g=richness_Traps_h_Arachnida$Mountain1Traps, p.adjust.method = "single-step")
##Comparison of each group against. 
axis(1, at=NULL, label=F, tick=TRUE, lwd=1.3)
text(x=c(1,2,3,4), y=(0.3), labels=c("a","b","a","a"), cex=2)
text(x=4.3, y=0.3, labels="*", cex=3)


read.table("../genetic/Data_out/Arachnida/Arachnida3P/richness_Traps0.03_Arachnida.txt",header=TRUE)->richness_Traps0.03

plot(richness_Traps0.03$Mountain1Traps,richness_Traps0.03$sample_richness_Traps0.03, xaxt="n", ylim=c(0,8), cex=1.4, cex.axis=2.3, lwd=2.5)
kruskal.test(sample_richness_Traps0.03 ~ Mountain1Traps, data = richness_Traps0.03)
PMCMRplus::kwAllPairsConoverTest(x=richness_Traps0.03$sample_richness_Traps0.03, g=richness_Traps0.03$Mountain1Traps, p.adjust.method = "single-step")
##Comparison of each group against. 
axis(1, at=NULL, label=F, tick=TRUE, lwd=1.3)
#text(x=c(1,2,3,4.1), y=(49), labels=c("ab","a","ab","b"), cex=1.4)
#text(x=c(1,2,3,4), y=(2), labels=c("ab","a","ab","b"), cex=1.5)
text(x=4.3, y=0.3, labels="ns", cex=2)

read.table("../genetic/Data_out/Arachnida/Arachnida5P/richness_Traps0.05_Arachnida.txt",header=TRUE)->richness_Traps0.05_Arachnida

plot(richness_Traps0.05_Arachnida$Mountain1Traps,richness_Traps0.05_Arachnida$sample_richness_Traps0.05, xaxt="n", ylim=c(0,8), cex=1.4, cex.axis=2.3, lwd=2.5)
kruskal.test(sample_richness_Traps0.05 ~ Mountain1Traps, data = richness_Traps0.05_Arachnida)
#posthoc.kruskal.nemenyi.test(x=richness_Traps0.05$sample_richness_Traps0.05, g=richness_Traps0.05$Mountain1Traps, method="Bonferroni")
PMCMRplus::kwAllPairsConoverTest(x=richness_Traps0.05_Arachnida$sample_richness_Traps0.05, g=richness_Traps0.05_Arachnida$Mountain1Traps, p.adjust.method = "single-step")
##Comparison of each group against. 
axis(1, at=NULL, label=F, tick=TRUE, lwd=1.3)
#text(x=c(1,2,3,4.1), y=(49), labels=c("ab","a","ab","b"), cex=1.4)
#text(x=c(1,2,3,4), y=(2), labels=c("ab","a","ab","b"), cex=1.5)
text(x=4.3, y=0.3, labels="ns", cex=2)
mtext(c("Arachnida"), side = 4, col = "black", line = 2, cex = 2)

#


#**Hempitera**
read.table("../genetic/Data_out/Hemiptera/Hemiptera_Haplotypes/richness_Traps_h_Hemiptera.txt",header=TRUE)->richness_Traps_h_Hemiptera

plot(richness_Traps_h_Hemiptera$Mountain1Traps,richness_Traps_h_Hemiptera$sample_richness_Traps_h, xaxt="n", cex=1.4, cex.axis=2.3, lwd=2.5, ylim=c(0,20))
kruskal.test(sample_richness_Traps_h ~ Mountain1Traps, data = richness_Traps_h_Hemiptera)
PMCMRplus::kwAllPairsConoverTest(x=richness_Traps_h_Hemiptera$sample_richness_Traps_h, g=richness_Traps_h_Hemiptera$Mountain1Traps, p.adjust.method = "single-step")
##Comparison of each group against.
axis(1, at=NULL, label=F, tick=TRUE, lwd=1.3)
#text(x=c(1,2,3,4), y=(2), labels=c("ab","a","ab","b"), cex=1.5)
text(x=4.3, y=0.3, labels="ns", cex=2)


read.table("../genetic/Data_out/Hemiptera/Hemiptera3P/richness_Traps0.03_Hemiptera.txt",header=TRUE)->richness_Traps0.03

plot(richness_Traps0.03$Mountain1Traps,richness_Traps0.03$sample_richness_Traps0.03, xaxt="n", ylim=c(0,15), cex=1.4, cex.axis=2.3, lwd=2.5)
kruskal.test(sample_richness_Traps0.03 ~ Mountain1Traps, data = richness_Traps0.03)
PMCMRplus::kwAllPairsConoverTest(x=richness_Traps0.03$sample_richness_Traps0.03, g=richness_Traps0.03$Mountain1Traps, p.adjust.method = "single-step")
##Comparison of each group against. 
axis(1, at=NULL, label=F, tick=TRUE, lwd=1.3)
#text(x=c(1,2,3,4.1), y=(49), labels=c("ab","a","ab","b"), cex=1.4)
#text(x=c(1,2,3,4), y=(2), labels=c("ab","a","ab","b"), cex=1.5)
text(x=4.3, y=0.3, labels="ns", cex=2)


read.table("../genetic/Data_out/Hemiptera/Hemiptera5P/richness_Traps0.05_Hemiptera.txt",header=TRUE)->richness_Traps0.05

plot(richness_Traps0.05$Mountain1Traps,richness_Traps0.05$sample_richness_Traps0.05, xaxt="n", ylim=c(0,15), cex=1.4, cex.axis=2.3, lwd=2.5)
kruskal.test(sample_richness_Traps0.05 ~ Mountain1Traps, data = richness_Traps0.05)
PMCMRplus::kwAllPairsConoverTest(x=richness_Traps0.05$sample_richness_Traps0.05, g=richness_Traps0.05$Mountain1Traps, p.adjust.method = "single-step")
##Comparison of each group against. 
axis(1, at=NULL, label=F, tick=TRUE, lwd=1.3)
#text(x=c(1,2,3,4.1), y=(49), labels=c("ab","a","ab","b"), cex=1.4)
#text(x=c(1,2,3,4), y=(2), labels=c("ab","a","ab","b"), cex=1.5)
text(x=4.3, y=0.3, labels="ns", cex=2)
mtext(c("Hemiptera"), side = 4, col = "black", line = 2, cex=2)

#

#**Hymenoptera**

read.table("../genetic/Data_out/Hymenoptera/Hymenoptera_Haplotypes/richness_Traps_h_Hymenoptera.txt",header=TRUE)->richness_Traps_h_Hymenoptera

plot(richness_Traps_h_Hymenoptera$Mountain1Traps,richness_Traps_h_Hymenoptera$sample_richness_Traps_h, xaxt="n", cex=1.4, cex.axis=2.3, lwd=2.5, ylim=c(0,15))
kruskal.test(sample_richness_Traps_h ~ Mountain1Traps, data = richness_Traps_h_Hymenoptera)
PMCMRplus::kwAllPairsConoverTest(x=richness_Traps_h_Hymenoptera$sample_richness_Traps_h, g=richness_Traps_h_Hymenoptera$Mountain1Traps, p.adjust.method = "single-step")
##Comparison of each group against. 
axis(1, at=NULL, label=F, tick=TRUE, lwd=1.3)
#text(x=c(1,2,3,4), y=(2), labels=c("ab","a","ab","b"), cex=1.5)
text(x=4.3, y=0.3, labels="ns", cex=2)


read.table("../genetic/Data_out/Hymenoptera/Hymenoptera3P/richness_Traps0.03_Hymenoptera.txt",header=TRUE)->richness_Traps0.03

plot(richness_Traps0.03$Mountain1Traps,richness_Traps0.03$sample_richness_Traps0.03, xaxt="n", ylim=c(0,10), cex=1.4, cex.axis=2.3, lwd=2.5)
kruskal.test(sample_richness_Traps0.03 ~ Mountain1Traps, data = richness_Traps0.03)
PMCMRplus::kwAllPairsConoverTest(x=richness_Traps0.03$sample_richness_Traps0.03, g=richness_Traps0.03$Mountain1Traps, p.adjust.method = "single-step")
##Comparison of each group against. 
axis(1, at=NULL, label=F, tick=TRUE, lwd=1.3)
#text(x=c(1,2,3,4.1), y=(49), labels=c("ab","a","ab","b"), cex=1.4)
#text(x=c(1,2,3,4), y=(2), labels=c("ab","a","ab","b"), cex=1.5)
text(x=4.3, y=0.3, labels="ns", cex=2)

read.table("../genetic/Data_out/Hymenoptera/Hymenoptera5P/richness_Traps0.05_Hymenoptera.txt",header=TRUE)->richness_Traps0.05

plot(richness_Traps0.05$Mountain1Traps,richness_Traps0.05$sample_richness_Traps0.05, xaxt="n", ylim=c(0,10), cex=1.4, cex.axis=2.3, lwd=2.5)
kruskal.test(sample_richness_Traps0.05 ~ Mountain1Traps, data = richness_Traps0.05)
PMCMRplus::kwAllPairsConoverTest(x=richness_Traps0.05$sample_richness_Traps0.05, g=richness_Traps0.05$Mountain1Traps, p.adjust.method = "single-step")
##Comparison of each group against. 
axis(1, at=NULL, label=F, tick=TRUE, lwd=1.3)
#text(x=c(1,2,3,4.1), y=(49), labels=c("ab","a","ab","b"), cex=1.4)
#text(x=c(1,2,3,4), y=(2), labels=c("ab","a","ab","b"), cex=1.5)
text(x=4.3, y=0.3, labels="ns", cex=2)
mtext(c("Hymenoptera"), side = 4, col = "black", line = 2, cex = 2)
#

#**Coleoptera**

read.table("../genetic/Data_out/Coleoptera/Coleoptera_Haplotypes/richness_Traps_h_Coleoptera.txt",header=TRUE)->richness_Traps_h_Coleoptera

plot(richness_Traps_h_Coleoptera$Mountain1Traps,richness_Traps_h_Coleoptera$sample_richness_Traps_h, xaxt="n", cex=1.4, cex.axis=2.3, lwd=2.5, ylim=c(0,6))
kruskal.test(sample_richness_Traps_h ~ Mountain1Traps, data = richness_Traps_h_Coleoptera)
PMCMRplus::kwAllPairsConoverTest(x=richness_Traps_h_Coleoptera$sample_richness_Traps_h, g=richness_Traps_h_Coleoptera$Mountain1Traps, p.adjust.method = "single-step")
##Comparison of each group against. 
axis(1, at=NULL, label=F, tick=TRUE, lwd=1.3)
#text(x=c(1,2,3,4), y=(2), labels=c("ab","a","ab","b"), cex=1.5)
text(x=c(1,2,3,4), -0.9, labels=c("S","P","B","W"), cex=2, xpd=NA)
text(x=4.3, y=0.3, labels="ns", cex=2)

read.table("../genetic/Data_out/Coleoptera/Coleoptera3P/richness_Traps0.03_Coleoptera.txt",header=TRUE)->richness_Traps0.03


plot(richness_Traps0.03$Mountain1Traps,richness_Traps0.03$sample_richness_Traps0.03, xaxt="n", ylim=c(0,5), cex=1.4, cex.axis=2.3, lwd=2.5)
kruskal.test(sample_richness_Traps0.03 ~ Mountain1Traps, data = richness_Traps0.03)
PMCMRplus::kwAllPairsConoverTest(x=richness_Traps0.03$sample_richness_Traps0.03, g=richness_Traps0.03$Mountain1Traps, p.adjust.method = "single-step")
##Comparison of each group against. 
axis(1, at=NULL, label=F, tick=TRUE, lwd=1.3)
#text(x=c(1,2,3,4.1), y=(49), labels=c("ab","a","ab","b"), cex=1.4)
#text(x=c(1,2,3,4), y=(2), labels=c("ab","a","ab","b"), cex=1.5)
text(x=c(1,2,3,4), -0.9, labels=c("S","P","B","W"), cex=2, xpd=NA)
text(x=4.3, y=0.3, labels="ns", cex=2)

read.table("../genetic/Data_out/Coleoptera/Coleoptera5P/richness_Traps0.05_Coleoptera.txt",header=TRUE)->richness_Traps0.05

plot(richness_Traps0.05$Mountain1Traps,richness_Traps0.05$sample_richness_Traps0.05, xaxt="n", ylim=c(0,5), cex=1.4, cex.axis=2.3, lwd=2.5)
kruskal.test(sample_richness_Traps0.05 ~ Mountain1Traps, data = richness_Traps0.05)
PMCMRplus::kwAllPairsConoverTest(x=richness_Traps0.05$sample_richness_Traps0.05, g=richness_Traps0.05$Mountain1Traps, p.adjust.method = "single-step")
##Comparison of each group against. 
axis(1, at=NULL, label=F, tick=TRUE, lwd=1.3)
#text(x=c(1,2,3,4.1), y=(49), labels=c("ab","a","ab","b"), cex=1.4)
#text(x=c(1,2,3,4), y=(2), labels=c("ab","a","ab","b"), cex=1.5)
text(x=4.3, y=0.3, labels="ns", cex=2)
text(x=c(1,2,3,4), -0.9, labels=c("S","P","B","W"), cex=2, xpd=NA)
mtext(c("Coleoptera"), side = 4, col = "black", line = 2, cex = 2)

##Legend
mtext("Type of Traps", side=1, outer=TRUE, line=2, cex=2.5)
mtext("R i c h n e s s", side=2, outer=TRUE, line=0.7, cex=2.5) 

dev.off()
#

#**plot Global Richness Experiments**
png(filename="../figures/Figure2_Richness_Experiment.png", width=620 , height=430, units="px") # set size of the file to plot 
par(mfrow=c(2,3), mar = c(3, 2.5, 2, 2), omi=c(0.5, 0.5, 0.5, 0.5)) #the number of rows and columns the figure would have

##**Creating plot**
##**Community diversity and composition at the haplotype, 3% and 5% lineages for each of the eight taxonomic orders studied**

#####This script gets all scripts of diversity using 6 arthropods order at multi-hierarchical levels 
#####We used 6 groups: Arachnda, Coleoptera, Collembola, Diptera, Hemiptera, and Hymenoptera.

#**plot Global Richness by Sites**
#**Diptera**

read.table("../genetic/Data_out/Diptera/Diptera_Haplotypes/richness_Groups_h_Diptera.txt",header=TRUE)->richness_Groups_h_Diptera

plot(richness_Groups_h_Diptera$Mountain2Groups,richness_Groups_h_Diptera$sample_richness_Groups_h, xaxt="n", cex=1.4, cex.axis=2.3, lwd=2.5, ylim=c(0,60), col= c("#153a7b", "#eaa22f", "#81b9d0"))
kruskal.test(sample_richness_Groups_h ~ Mountain2Groups, data = richness_Groups_h_Diptera)
PMCMRplus::kwAllPairsConoverTest(x=richness_Groups_h_Diptera$sample_richness_Groups_h, g=richness_Groups_h_Diptera$Mountain2Groups) 
##Comparison of each group against. 
axis(1, at=NULL, label=F, tick=TRUE, lwd=1.3)
#text(x=c(1,2,3,4), y=(2), labels=c("ab","a","ab","b"), cex=2)
text(x=3.3, y=1, labels="ns", cex=2)
mtext(c("haplotypes"), side = 3, col = "black", line = 1, cex=1.5)

read.table("../genetic/Data_out/Diptera/Diptera3P/richness_Groups0.03_Diptera.txt",header=TRUE)->richness_Groups0.03

plot(richness_Groups0.03$Mountain2Groups,richness_Groups0.03$sample_richness_Groups0.03, xaxt="n", cex=1.4, cex.axis=2.3, lwd=2.5, ylim=c(0,30), col= c("#153a7b", "#eaa22f", "#81b9d0"))
kruskal.test(sample_richness_Groups0.03 ~ Mountain2Groups, data = richness_Groups0.03)
PMCMRplus::kwAllPairsConoverTest(x=richness_Groups0.03$sample_richness_Groups0.03, g=richness_Groups0.03$Mountain2Groups, p.adjust.method = "single-step")
##Comparison of each group against. 
axis(1, at=NULL, label=F, tick=TRUE, lwd=1.3)
text(x=c(1,2,3), y=1, labels=c("a","ab","b"), cex=2)
text(x=3.3, y=1, labels="*", cex=3.5)
mtext(c("3% lineages"), side = 3, col = "black", line = 1, cex=1.5)

read.table("../genetic/Data_out/Diptera/Diptera5P/richness_Groups0.05_Diptera.txt",header=TRUE)->richness_Groups0.05

plot(richness_Groups0.05$Mountain2Groups,richness_Groups0.05$sample_richness_Groups0.05, xaxt="n", cex=1.4, cex.axis=2.3, lwd=2.5, ylim=c(0,30), col= c("#153a7b", "#eaa22f", "#81b9d0"))
kruskal.test(sample_richness_Groups0.05 ~ Mountain2Groups, data = richness_Groups0.05)
PMCMRplus::kwAllPairsConoverTest(x=richness_Groups0.05$sample_richness_Groups0.05, g=richness_Groups0.05$Mountain2Groups, p.adjust.method = "single-step")
axis(1, at=NULL, label=F, tick=TRUE, lwd=1.3)
text(x=c(1,2,3), y=1, labels=c("a","ab","b"), cex=2)
text(x=3.3, y=1, labels="*", cex=3.5)
mtext(c("5% lineages"), side = 3, col = "black", line = 1, cex = 1.5)
mtext(c("Diptera"), side = 4, col = "black", line = 2, cex=1.5)

#**Hymenoptera**

read.table("../genetic/Data_out/Hymenoptera/Hymenoptera_Haplotypes/richness_Groups_h_Hymenoptera.txt",header=TRUE)->richness_Groups_h_Hymenoptera

plot(richness_Groups_h_Hymenoptera$Mountain2Groups,richness_Groups_h_Hymenoptera$sample_richness_Groups_h, xaxt="n", cex=1.4, cex.axis=2.3, lwd=2.5, ylim=c(0,20), col= c("#153a7b", "#eaa22f", "#81b9d0"))
kruskal.test(sample_richness_Groups_h ~ Mountain2Groups, data = richness_Groups_h_Hymenoptera)
PMCMRplus::kwAllPairsConoverTest(x=richness_Groups_h_Hymenoptera$sample_richness_Groups_h, g=richness_Groups_h_Hymenoptera$Mountain2Groups, p.adjust.method = "single-step")
axis(1, at=NULL, label=F, tick=TRUE, lwd=1.3)
text(x=3.3, y=0.5, labels="ns", cex=2)
text(x=c(1,2,3), -2.2, labels=c("A","B","AB"), cex=1.6, xpd=NA)
#mtext(c("haplotypes"), side = 3, col = "black", line = 1, cex=2)


read.table("../genetic/Data_out/Hymenoptera/Hymenoptera3P/richness_Groups0.03_Hymenoptera.txt",header=TRUE)->richness_Groups0.03

plot(richness_Groups0.03$Mountain2Groups,richness_Groups0.03$sample_richness_Groups0.03, xaxt="n", cex=1.4, cex.axis=2.3, lwd=2.5, ylim=c(0,15), col= c("#153a7b", "#eaa22f", "#81b9d0"))
kruskal.test(sample_richness_Groups0.03 ~ Mountain2Groups, data = richness_Groups0.03)
PMCMRplus::kwAllPairsConoverTest(x=richness_Groups0.03$sample_richness_Groups0.03, g=richness_Groups0.03$Mountain2Groups, p.adjust.method = "single-step")
##Comparison of each group against. 
axis(1, at=NULL, label=F, tick=TRUE, lwd=1.3)
#text(x=c(1,2,3,4,5), y=1.4, labels=c("ab","ab","ab","a","b"), cex=1.5)
text(x=3.3, y=0.5, labels="ns", cex=2)
text(x=c(1,2,3), -2, labels=c("A","B","AB"), cex=1.6, xpd=NA)
#mtext(c("lineages 3%"), side = 3, col = "black", line = 1, cex = 2)

read.table("../genetic/Data_out/Hymenoptera/Hymenoptera5P/richness_Groups0.05_Hymenoptera.txt",header=TRUE)->richness_Groups0.05

plot(richness_Groups0.05$Mountain2Groups,richness_Groups0.05$sample_richness_Groups0.05, xaxt="n", cex=1.4, cex.axis=2.3, lwd=2.5, ylim=c(0,15),col= c("#153a7b", "#eaa22f", "#81b9d0"))
kruskal.test(sample_richness_Groups0.05 ~ Mountain2Groups, data = richness_Groups0.05)
PMCMRplus::kwAllPairsConoverTest(x=richness_Groups0.05$sample_richness_Groups0.05, g=richness_Groups0.05$Mountain2Groups, p.adjust.method = "single-step")
axis(1, at=NULL, label=F, tick=TRUE, lwd=1.3)
text(x=3.3, y=0.5, labels="ns", cex=2)
#mtext(c("lineages 5%"), side = 3, col = "black", line = 1, cex = 2)
text(x=c(1,2,3), -2, labels=c("A","AB","AB"), cex=1.6, xpd=NA)
mtext(c("Hymenoptera"), side = 4, col = "black", line = 2, cex = 1.5)

##Legend
mtext("Experiment Samples", side=1, outer=TRUE, line=2, cex=1.8)
mtext("R i c h n e s s", side=2, outer=TRUE, line=0.7, cex=1.8) 

dev.off()
#


#**Non-Metric Multidimensional scaling (NMDS) ordinations of community similarity**

png(filename="../figures/Figure3_Non_MetricMultidimensionalScaling.png", width=620, height=460, units="px") # set size of the file to plot 
par(mfrow=c(2,3), mar = c(1.8, 1.8, 1.8, 1.8), omi=c(0.3, 0.3, 0.3, 0.3)) #the number of rows and columns the figure would have

#**Creating plot**
#**Diptera**
read.table ("../genetic/Data_out/Diptera/Diptera_Haplotypes/community_Diptera_Groups_h.txt")->community_Diptera_Groups_h
read.table("../genetic/Data_out/Diptera/Diptera_Haplotypes/general_sample_Mountain2Groups_h.txt",header=TRUE)->general_sample_Mountain2Groups_h

beta.multi(community_Diptera_Groups_h, index.family="sorensen")
##turnover by pairs, nmds, anosim
beta.pair(community_Diptera_Groups_h, index.family="sorensen")->beta.pair  ##betadiversity by pair of communities using sorensen on the precense/absence data, with estimation of turnover and nestedness datamatrixes simultaneously
metaMDS (beta.pair$beta.sim)->MDSbetasim_Gh ##NMDS
plot (MDSbetasim_Gh, xlim=c(-0.5, 0.5), ylim=c(-0.5, 0.5), cex=2, cex.lab=1, cex.axis=2.3, lwd=4.5)
with(general_sample_Mountain2Groups_h,ordispider(MDSbetasim_Gh, Mountain2, cex.lab=1, col= c("#153a7b", "#eaa22f", "#81b9d0"), lwd=4.5))
##I put letter "r" in cursive and r2 value
mylabel = bquote(italic(r)^2)
text(x=0.10, y=-0.45, labels = mylabel, cex=1.5)
text(x=0.35, y=-0.45, labels="=0.129 ns", cex=1.5)
mtext(c("Haplotypes"), side = 3, col = "black", line = 1, cex = 1.7)


read.table ("../genetic/Data_out/Diptera/Diptera3P/community_Diptera_Groups0.03.txt")->community_Diptera_Groups0.03
read.table("../genetic/Data_out/Diptera/Diptera3P/general_sample_Mountain2Groups0.03.txt",header=TRUE)->general_sample_Mountain2Groups0.03

beta.multi(community_Diptera_Groups0.03, index.family="sorensen")
##turnover by pairs, nmds, anosim
beta.pair(community_Diptera_Groups0.03, index.family="sorensen")->beta.pair  ##betadiversity by pair of communities using sorensen on the precense/absence data, with estimation of turnover and nestedness datamatrixes simultaneously
metaMDS (beta.pair$beta.sim)->MDSbetasim_G0.03 ##NMDS
plot (MDSbetasim_G0.03, xlim=c(-0.5, 0.5), ylim=c(-0.5, 0.5), cex=2, cex.lab=1, cex.axis=2.3, lwd=4.5)
with(general_sample_Mountain2Groups0.03,ordispider(MDSbetasim_G0.03, Mountain2, cex.lab=1, col= c("#153a7b", "#eaa22f", "#81b9d0"), lwd=4.5))
##I put letter "r" in cursive and r2 value
mylabel = bquote(italic(r)^2)
text(x=0.10, y=-0.45, labels = mylabel, cex=1.5)
text(x=0.35, y=-0.45, labels="=0.24 ns", cex=1.5)
mtext(c("lineages 3%"), side = 3, col = "black", line = 1, cex = 1.7)

read.table ("../genetic/Data_out/Diptera/Diptera5P/community_Diptera_Groups0.05.txt")->community_Diptera_Groups0.05
read.table("../genetic/Data_out/Diptera/Diptera5P/general_sample_Mountain2Groups0.05.txt",header=TRUE)->general_sample_Mountain2Groups0.05

beta.multi(community_Diptera_Groups0.05, index.family="sorensen")
##turnover by pairs, nmds, anosim
beta.pair(community_Diptera_Groups0.05, index.family="sorensen")->beta.pair  ##betadiversity by pair of communities using sorensen on the precense/absence data, with estimation of turnover and nestedness datamatrixes simultaneously
metaMDS (beta.pair$beta.sim)->MDSbetasim_G0.05 ##NMDS
plot (MDSbetasim_G0.05, xlim=c(-0.5, 0.5), ylim=c(-0.5, 0.5), cex=2, cex.lab=1.4, cex.axis=2.3, lwd=4.5)
with(general_sample_Mountain2Groups0.05,ordispider(MDSbetasim_G0.05, Mountain2, cex.lab=1, col= c("#153a7b", "#eaa22f", "#81b9d0"),lwd=4.5))
##I put letter "r" in cursive and r2 value
mylabel = bquote(italic(r)^2)
text(x=0.10, y=-0.45, labels = mylabel, cex=1.5)
text(x=0.35, y=-0.45, labels="=0.29 ns", cex=1.5)
mtext(c("lineages 5%"), side = 3, col = "black", line = 1, cex = 1.7)
mtext(c("Diptera"), side = 4, col = "black", line = 2, cex = 1.7)
#


#**Hymenoptera**

read.table ("../genetic/Data_out/Hymenoptera/Hymenoptera_Haplotypes/community_Hymenoptera_Groups_h.txt")->community_Hymenoptera_Groups_h
read.table("../genetic/Data_out/Hymenoptera/Hymenoptera_Haplotypes/general_sample_Mountain2Groups_h.txt",header=TRUE)->general_sample_Mountain2Groups_h

beta.multi(community_Hymenoptera_Groups_h, index.family="sorensen")
##turnover by pairs, nmds, anosim
beta.pair(community_Hymenoptera_Groups_h, index.family="sorensen")->beta.pair  ##betadiversity by pair of communities using sorensen on the precense/absence data, with estimation of turnover and nestedness datamatrixes simultaneously
metaMDS (beta.pair$beta.sim)->MDSbetasim_Gh ##NMDS
plot (MDSbetasim_Gh, xlim=c(-0.5, 0.5), ylim=c(-0.5, 0.5), cex=2, cex.lab=1, cex.axis=2.3, lwd=4.5)
with(general_sample_Mountain2Groups_h,ordispider(MDSbetasim_Gh, Mountain2, cex.lab=1, col= c("#153a7b", "#eaa22f", "#81b9d0"), lwd=4.5))
##I put letter "r" in cursive and r2 value
mylabel = bquote(italic(r)^2)
text(x=0.10, y=-0.45, labels = mylabel, cex=1.5)
text(x=0.35, y=-0.45, labels="=0.008 ns", cex=1.5)


read.table ("../genetic/Data_out/Hymenoptera/Hymenoptera3P/community_Hymenoptera_Groups0.03.txt")->community_Hymenoptera_Groups0.03
read.table("../genetic/Data_out/Hymenoptera/Hymenoptera3P/general_sample_Mountain2Groups0.03.txt",header=TRUE)->general_sample_Mountain2Groups0.03

beta.multi(community_Hymenoptera_Groups0.03, index.family="sorensen")
##turnover by pairs, nmds, anosim
beta.pair(community_Hymenoptera_Groups0.03, index.family="sorensen")->beta.pair  ##betadiversity by pair of communities using sorensen on the precense/absence data, with estimation of turnover and nestedness datamatrixes simultaneously
metaMDS (beta.pair$beta.sim)->MDSbetasim_G0.03 ##NMDS
plot (MDSbetasim_G0.03, xlim=c(-0.5, 0.5), ylim=c(-0.5, 0.5), cex=2, cex.lab=1, cex.axis=2.3, lwd=4.5)
with(general_sample_Mountain2Groups0.03,ordispider(MDSbetasim_G0.03, Mountain2, cex.lab=1, col= c("#153a7b", "#eaa22f", "#81b9d0"), lwd=4.5))
##I put letter "r" in cursive and r2 value
mylabel = bquote(italic(r)^2)
text(x=0.10, y=-0.45, labels = mylabel, cex=1.5)
text(x=0.35, y=-0.45, labels="=0.026 ns", cex=1.5)

read.table ("../genetic/Data_out/Hymenoptera/Hymenoptera5P/community_Hymenoptera_Groups0.05.txt")->community_Hymenoptera_Groups0.05
read.table("../genetic/Data_out/Hymenoptera/Hymenoptera5P/general_sample_Mountain2Groups0.05.txt",header=TRUE)->general_sample_Mountain2Groups0.05

beta.multi(community_Hymenoptera_Groups0.05, index.family="sorensen")
##turnover by pairs, nmds, anosim
beta.pair(community_Hymenoptera_Groups0.05, index.family="sorensen")->beta.pair  ##betadiversity by pair of communities using sorensen on the precense/absence data, with estimation of turnover and nestedness datamatrixes simultaneously
metaMDS (beta.pair$beta.sim)->MDSbetasim_G0.05 ##NMDS
plot (MDSbetasim_G0.05, xlim=c(-0.5, 0.5), ylim=c(-0.5, 0.5), cex=2, cex.lab=1, cex.axis=2.3, lwd=4.5)
with(general_sample_Mountain2Groups0.05,ordispider(MDSbetasim_G0.05, Mountain2, cex.lab=1, col= c("#153a7b", "#eaa22f", "#81b9d0"), lwd=4.5))
##I put letter "r" in cursive and r2 value
mylabel = bquote(italic(r)^2)
text(x=0.10, y=-0.45, labels = mylabel, cex=1.5)
text(x=0.35, y=-0.45, labels="=0.027 ns", cex=1.5)
mtext(c("Hymenoptera"), side = 4, col = "black", line = 2, cex = 2)

## Legend
#legend(0.8, 0.4, xpd=NA,
#       legend = c("DNA from large", "DNA from small", "pool", "Tlacotepec"),
#       col=c("#153a7b", "#eaa22f", "#97518b", "#81b9d0"),
#       bty="n", text.font=2, lty=1, lwd= 8, cex = 2.5)

dev.off()

#



#**Accumulation Curves**
png(filename="../figures/FigureS3_AccumulationCurves.png", width=920 , height=1536, units="px") # set size of the file to plot 
par(mfrow=c(8,3), mar = c(2, 2, 2, 2), omi=c(0.5, 0.5, 0.5, 0.5)) #the number of rows and columns the figure would have

#**Creating plot**
#**Diptera**
read.table("../genetic/Data_out/Diptera/Diptera_Haplotypes/community_Diptera_Site_h.txt")->community_Diptera_Site_h

specaccum(community_Diptera_Site_h,"random", permutations=1000)->cum_Site_h
plot(cum_Site_h, cex=1.4, cex.lab=1.4, cex.axis=2.3, lwd=3, xlim=c(0,42))
specpool(community_Diptera_Site_h)->specpool_Site_h
specpool_Site_h$Species/specpool_Site_h$chao*100
text(x=39, y=3, labels="87.66%", cex=2)
mtext(c("haplotypes"), side = 3, col = "black", line = 1, cex = 2)

read.table("../genetic/Data_out/Diptera/Diptera3P/community_Diptera_Site0.03.txt")->community_Diptera_Site0.03

specaccum(community_Diptera_Site0.03,"random", permutations=1000)->cum_Site0.03
plot(cum_Site0.03, cex=1.4, cex.lab=1.4, cex.axis=2.3, lwd=3, ylim=c(0,180), xlim=c(0,42))
specpool(community_Diptera_Site0.03)->specpool_Site0.03
specpool_Site0.03$Species/specpool_Site0.03$chao*100
text(x=39, y=3, labels="80.66%", cex=2)
mtext(c("3% lineages"), side = 3, col = "black", line = 1, cex = 2)

read.table("../genetic/Data_out/Diptera/Diptera5P/community_Diptera_Site0.05.txt")->community_Diptera_Site0.05

specaccum(community_Diptera_Site0.05,"random", permutations=1000)->cum_Site0.05
plot(cum_Site0.05, cex=1.4, cex.lab=1.4, cex.axis=2.3, lwd=3, ylim=c(0,180), xlim=c(0,42))
specpool(community_Diptera_Site0.05)->specpool_Site0.05
specpool_Site0.05$Species/specpool_Site0.05$chao*100
text(x=39, y=3, labels="82.84%", cex=2)
mtext(c("5% lineages"), side = 3, col = "black", line = 1, cex = 2)
mtext(c("Diptera"), side = 4, col = "black", line = 2, cex = 2)
#

#**Collembola**
read.table("../genetic/Data_out/Collembola/Collembola_Haplotypes/community_Collembola_Site_h.txt")->community_Collembola_Site_h

specaccum(community_Collembola_Site_h,"random", permutations=1000)->cum_Site_h
plot(cum_Site_h, cex=1.4, cex.lab=1.4, cex.axis=2.3, lwd=3, xlim=c(0,42))
specpool(community_Collembola_Site_h)->specpool_Site_h
specpool_Site_h$Species/specpool_Site_h$chao*100
text(x=39, y=3, labels="91.05%", cex=2)

read.table("../genetic/Data_out/Collembola/Collembola3P/community_Collembola_Site0.03.txt")->community_Collembola_Site0.03

specaccum(community_Collembola_Site0.03,"random", permutations=1000)->cum_Site0.03
plot(cum_Site0.03, cex=1.4, cex.lab=1.4, cex.axis=2.3, lwd=3, ylim=c(0,100), xlim=c(0,42))
specpool(community_Collembola_Site0.03)->specpool_Site0.03
specpool_Site0.03$Species/specpool_Site0.03$chao*100
text(x=39, y=1.5, labels="85.00%", cex=2)

read.table("../genetic/Data_out/Collembola/Collembola5P/community_Collembola_Site0.05.txt")->community_Collembola_Site0.05

specaccum(community_Collembola_Site0.05,"random", permutations=1000)->cum_Site0.05
plot(cum_Site0.05, cex=1.4, cex.lab=1.4, cex.axis=2.3, lwd=3, ylim=c(0,100), xlim=c(0,42))
specpool(community_Collembola_Site0.05)->specpool_Site0.05
specpool_Site0.05$Species/specpool_Site0.05$chao*100
text(x=39, y=1.5, labels="76.90%", cex=2)
mtext(c("Collembola"), side = 4, col = "black", line = 2, cex = 2)
#

#**Arachnida**
read.table("../genetic/Data_out/Arachnida/Arachnida_Haplotypes/community_Arachnida_Site_h.txt")->community_Arachnida_Site_h

specaccum(community_Arachnida_Site_h,"random", permutations=1000)->cum_Site_h
plot(cum_Site_h, cex=1.4, cex.lab=1.4, cex.axis=2.3, lwd=3, xlim=c(0,42))
specpool(community_Arachnida_Site_h)->specpool_Site_h
specpool_Site_h$Species/specpool_Site_h$chao*100
text(x=39, y=3, labels="66.99%", cex=2)

read.table("../genetic/Data_out/Arachnida/Arachnida3P/community_Arachnida_Site0.03.txt")->community_Arachnida_Site0.03

specaccum(community_Arachnida_Site0.03,"random", permutations=1000)->cum_Site0.03
plot(cum_Site0.03, cex=1.4, cex.lab=1.4, cex.axis=2.3, lwd=3, ylim=c(0,100), xlim=c(0,42))
specpool(community_Arachnida_Site0.03)->specpool_Site0.03
specpool_Site0.03$Species/specpool_Site0.03$chao*100
text(x=39, y=1.5, labels="60.59%", cex=2)

read.table("../genetic/Data_out/Arachnida/Arachnida5P/community_Arachnida_Site0.05.txt")->community_Arachnida_Site0.05

specaccum(community_Arachnida_Site0.05,"random", permutations=1000)->cum_Site0.05
plot(cum_Site0.05, cex=1.4, cex.lab=1.4, cex.axis=2.3, lwd=3, ylim=c(0,100), xlim=c(0,42))
specpool(community_Arachnida_Site0.05)->specpool_Site0.05
specpool_Site0.05$Species/specpool_Site0.05$chao*100
text(x=39, y=1.5, labels="46.46%", cex=2)
mtext(c("Arachnida"), side = 4, col = "black", line = 2, cex = 2)
#

#**Hempitera**
read.table("../genetic/Data_out/Hemiptera/Hemiptera_Haplotypes/community_Hemiptera_Site_h.txt")->community_Hemiptera_Site_h

specaccum(community_Hemiptera_Site_h,"random", permutations=1000)->cum_Site_h
plot(cum_Site_h, cex=1.4, cex.lab=1.4, cex.axis=2.3, lwd=3, xlim=c(0,42))
specpool(community_Hemiptera_Site_h)->specpool_Site_h
specpool_Site_h$Species/specpool_Site_h$chao*100
text(x=39, y=3, labels="58.71%", cex=2)

read.table("../genetic/Data_out/Hemiptera/Hemiptera3P/community_Hemiptera_Site0.03.txt")->community_Hemiptera_Site0.03

specaccum(community_Hemiptera_Site0.03,"random", permutations=1000)->cum_Site0.03
plot(cum_Site0.03, cex=1.4, cex.lab=1.4, cex.axis=2.3, lwd=3, ylim=c(0,80), xlim=c(0,42))
specpool(community_Hemiptera_Site0.03)->specpool_Site0.03
specpool_Site0.03$Species/specpool_Site0.03$chao*100
text(x=39, y=1, labels="60.00%", cex=2)

read.table("../genetic/Data_out/Hemiptera/Hemiptera5P/community_Hemiptera_Site0.05.txt")->community_Hemiptera_Site0.05

specaccum(community_Hemiptera_Site0.05,"random", permutations=1000)->cum_Site0.05
plot(cum_Site0.05, cex=1.4, cex.lab=1.4, cex.axis=2.3, lwd=3, ylim=c(0,80), xlim=c(0,42))
specpool(community_Hemiptera_Site0.05)->specpool_Site0.05
specpool_Site0.05$Species/specpool_Site0.05$chao*100
text(x=39, y=1, labels="60.77%", cex=2)
mtext(c("Hemiptera"), side = 4, col = "black", line = 2, cex = 2)
#

#**Hymenoptera**
read.table("../genetic/Data_out/Hymenoptera/Hymenoptera_Haplotypes/community_Hymenoptera_Site_h.txt")->community_Hymenoptera_Site_h

specaccum(community_Hymenoptera_Site_h,"random", permutations=1000)->cum_Site_h
plot(cum_Site_h, ylim=c(0,135), cex=1.4, cex.lab=1.4, cex.axis=2.3, lwd=3, xlim=c(0,42))
specpool(community_Hymenoptera_Site_h)->specpool_Site_h
specpool_Site_h$Species/specpool_Site_h$chao*100
text(x=39, y=2, labels="64.83%", cex=2)

read.table("../genetic/Data_out/Hymenoptera/Hymenoptera3P/community_Hymenoptera_Site0.03.txt")->community_Hymenoptera_Site0.03

specaccum(community_Hymenoptera_Site0.03,"random", permutations=1000)->cum_Site0.03
plot(cum_Site0.03, cex=1.4, cex.lab=1.4, cex.axis=2.3, lwd=3, ylim=c(0,100), xlim=c(0,42))
specpool(community_Hymenoptera_Site0.03)->specpool_Site0.03
specpool_Site0.03$Species/specpool_Site0.03$chao*100
text(x=39, y=1.5, labels="66.69%", cex=2)

read.table("../genetic/Data_out/Hymenoptera/Hymenoptera5P/community_Hymenoptera_Site0.05.txt")->community_Hymenoptera_Site0.05

specaccum(community_Hymenoptera_Site0.05,"random", permutations=1000)->cum_Site0.05
plot(cum_Site0.05, cex=1.4, cex.lab=1.4, cex.axis=2.3, lwd=3, ylim=c(0,100), xlim=c(0,42))
specpool(community_Hymenoptera_Site0.05)->specpool_Site0.05
specpool_Site0.05$Species/specpool_Site0.05$chao*100
text(x=39, y=1.5, labels="69.42%", cex=2)
mtext(c("Hymenoptera"), side = 4, col = "black", line = 2, cex = 2)
#

#**Coleoptera**
read.table("../genetic/Data_out/Coleoptera/Coleoptera_Haplotypes/community_Coleoptera_Site_h.txt")->community_Coleoptera_Site_h

specaccum(community_Coleoptera_Site_h,"random", permutations=1000)->cum_Site_h
plot(cum_Site_h, cex=1.4, cex.lab=1.4, cex.axis=2.3, lwd=3, xlim=c(0,42))
specpool(community_Coleoptera_Site_h)->specpool_Site_h
specpool_Site_h$Species/specpool_Site_h$chao*100
text(x=39, y=3, labels="55.32%", cex=2)

read.table("../genetic/Data_out/Coleoptera/Coleoptera3P/community_Coleoptera_Site0.03.txt")->community_Coleoptera_Site0.03

specaccum(community_Coleoptera_Site0.03,"random", permutations=1000)->cum_Site0.03
plot(cum_Site0.03, cex=1.4, cex.lab=1.4, cex.axis=2.3, lwd=3, ylim=c(0,100), xlim=c(0,42))
specpool(community_Coleoptera_Site0.03)->specpool_Site0.03
specpool_Site0.03$Species/specpool_Site0.03$chao*100
text(x=39, y=1.5, labels="49.22%", cex=2)

read.table("../genetic/Data_out/Coleoptera/Coleoptera5P/community_Coleoptera_Site0.05.txt")->community_Coleoptera_Site0.05

specaccum(community_Coleoptera_Site0.05,"random", permutations=1000)->cum_Site0.05
plot(cum_Site0.05, cex=1.4, cex.lab=1.4, cex.axis=2.3, lwd=3, ylim=c(0,100), xlim=c(0,42))
specpool(community_Coleoptera_Site0.05)->specpool_Site0.05
specpool_Site0.05$Species/specpool_Site0.05$chao*100
text(x=39, y=1.5, labels="49.71%", cex=2)
mtext(c("Coleoptera"), side = 4, col = "black", line = 2, cex = 2)
#

#**Myriapoda**
read.table("../genetic/Data_out/Myriapoda/Myriapoda_Haplotypes/community_Myriapoda_Site_h.txt")->community_Myriapoda_Site_h

specaccum(community_Myriapoda_Site_h,"random", permutations=1000)->cum_Site_h
plot(cum_Site_h, ylim=c(0,65), xlim=c(0,42), cex=1.4, cex.lab=1.4, cex.axis=2.3, lwd=3)
specpool(community_Myriapoda_Site_h)->specpool_Site_h
specpool_Site_h$Species/specpool_Site_h$chao*100
text(x=39, y=1.5, labels="26.58%", cex=2)

read.table("../genetic/Data_out/Myriapoda/Myriapoda3P/community_Myriapoda_Site0.03.txt")->community_Myriapoda_Site0.03

specaccum(community_Myriapoda_Site0.03,"random", permutations=1000)->cum_Site0.03
plot(cum_Site0.03, cex=1.4, cex.lab=1.4, cex.axis=2.3, lwd=3, ylim=c(0,20), xlim=c(0,42))
specpool(community_Myriapoda_Site0.03)->specpool_Site0.03
specpool_Site0.03$Species/specpool_Site0.03$chao*100
text(x=39, y=0.3, labels="77.67%", cex=2)

read.table("../genetic/Data_out/Myriapoda/Myriapoda5P/community_Myriapoda_Site0.05.txt")->community_Myriapoda_Site0.05

specaccum(community_Myriapoda_Site0.05,"random", permutations=1000)->cum_Site0.05
plot(cum_Site0.05, cex=1.4, cex.lab=1.4, cex.axis=2.3, lwd=3, ylim=c(0,20), xlim=c(0,42))
specpool(community_Myriapoda_Site0.05)->specpool_Site0.05
specpool_Site0.05$Species/specpool_Site0.05$chao*100
text(x=39, y=0.3, labels="76.16%", cex=2)
mtext(c("Myriapoda"), side = 4, col = "black", line = 2, cex = 2)
#

#**Lepidoptera**
read.table("../genetic/Data_out/Lepidoptera/Lepidoptera_Haplotypes/community_Lepidoptera_Site_h.txt")->community_Lepidoptera_Site_h

specaccum(community_Lepidoptera_Site_h,"random", permutations=1000)->cum_Site_h
plot(cum_Site_h, cex=1.4, cex.lab=1.4, cex.axis=2.3, lwd=3, xlim=c(0,42))
specpool(community_Lepidoptera_Site_h)->specpool_Site_h
specpool_Site_h$Species/specpool_Site_h$chao*100
text(x=39, y=2, labels="57.14%", cex=2)

read.table("../genetic/Data_out/Lepidoptera/Lepidoptera3P/community_Lepidoptera_Site0.03.txt")->community_Lepidoptera_Site0.03

specaccum(community_Lepidoptera_Site0.03,"random", permutations=1000)->cum_Site0.03
plot(cum_Site0.03, cex=1.4, cex.lab=1.4, cex.axis=2.3, lwd=3, ylim=c(0,40), xlim=c(0,42))
specpool(community_Lepidoptera_Site0.03)->specpool_Site0.03
specpool_Site0.03$Species/specpool_Site0.03$chao*100
text(x=39, y=0.5, labels="50.36%", cex=2)

read.table ("../genetic/Data_out/Lepidoptera/Lepidoptera5P/community_Lepidoptera_Site0.05.txt")->community_Lepidoptera_Site0.05

specaccum(community_Lepidoptera_Site0.05,"random", permutations=1000)->cum_Site0.05
plot(cum_Site0.05, cex=1.4, cex.lab=1.4, cex.axis=2.3, lwd=3, ylim=c(0,40), xlim=c(0,42))
specpool(community_Lepidoptera_Site0.05)->specpool_Site0.05
specpool_Site0.05$Species/specpool_Site0.05$chao*100
text(x=39, y=0.5, labels="36.52%", cex=2)
mtext(c("Lepidoptera"), side = 4, col = "black", line = 2, cex = 2)

mtext("Sampled sites (biodiversity soup)", side=1, outer=TRUE, line=1.8, cex=2.3)
mtext("Species accumulation", side=2.3, outer=TRUE, line=1, cex=2.3) 

dev.off()

#**END**
