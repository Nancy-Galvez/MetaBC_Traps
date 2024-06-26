#!/bin/bash

#SBATCH --mem=64000
#SBATCH -n 10


#Arachnida_138_33SOPAS_UNIQUES_ExpeRelabel420_20_26_1_alignment.fasta

vsearch --search_exact ../genetic/Data_in/LENGTH420CatExpe33SOPAS.Ptrim.trimmo.Maxee1.LABEL.sorted420NOPuntoComaBUE.MOD2.fas --db ../genetic/AlignedSeq/Arachnida_138_33SOPAS_UNIQUES_ExpeRelabel420_20_26_1_alignment.fasta --otutabout ../genetic/Data_out/ASVs_Table/allele_table_Arachnida_138_33SOPAS_UNIQUES_ExpeRelabel420_20_26_1_alignment.txt


#Coleoptera_71_33SOPAS_UNIQUES_ExpeRelabel420_20_26_1_alignment.fasta

vsearch --search_exact ../genetic/Data_in/LENGTH420CatExpe33SOPAS.Ptrim.trimmo.Maxee1.LABEL.sorted420NOPuntoComaBUE.MOD2.fas --db ../genetic/AlignedSeq/Coleoptera_71_33SOPAS_UNIQUES_ExpeRelabel420_20_26_1_alignment.fasta --otutabout ../genetic/Data_out/ASVs_Table/allele_table_Coleoptera_71_33SOPAS_UNIQUES_ExpeRelabel420_20_26_1_alignment.txt 


#Collembola_158_33SOPAS_UNIQUES_ExpeRelabel420_20_26_1_alignment.fasta

vsearch --search_exact ../genetic/Data_in/LENGTH420CatExpe33SOPAS.Ptrim.trimmo.Maxee1.LABEL.sorted420NOPuntoComaBUE.MOD2.fas --db ../genetic/AlignedSeq/Collembola_158_33SOPAS_UNIQUES_ExpeRelabel420_20_26_1_alignment.fasta --otutabout ../genetic/Data_out/ASVs_Table/allele_table_Collembola_158_33SOPAS_UNIQUES_ExpeRelabel420_20_26_1_alignment.txt


#Dipteras_583_33SOPAS_UNIQUES_ExpeRelabel420_20_26_1_alignment.fasta

vsearch --search_exact ../genetic/Data_in/LENGTH420CatExpe33SOPAS.Ptrim.trimmo.Maxee1.LABEL.sorted420NOPuntoComaBUE.MOD2.fas --db ../genetic/AlignedSeq/Dipteras_583_33SOPAS_UNIQUES_ExpeRelabel420_20_26_1_alignment.fasta --otutabout ../genetic/Data_out/ASVs_Table/allele_table_Dipteras_583_33SOPAS_UNIQUES_ExpeRelabel420_20_26_1_alignment.txt


#Hemiptera_150_33SOPAS_UNIQUES_ExpeRelabel420_20_26_1_alignment.fasta

vsearch --search_exact ../genetic/Data_in/LENGTH420CatExpe33SOPAS.Ptrim.trimmo.Maxee1.LABEL.sorted420NOPuntoComaBUE.MOD2.fas --db ../genetic/AlignedSeq/Hemiptera_150_33SOPAS_UNIQUES_ExpeRelabel420_20_26_1_alignment.fasta --otutabout ../genetic/Data_out/ASVs_Table/allele_table_Hemiptera_150_33SOPAS_UNIQUES_ExpeRelabel420_20_26_1_alignment.txt 


#Hymenoptera_161_33SOPAS_UNIQUES_ExpeRelabel420_20_26_1_alignment.fasta

vsearch --search_exact ../genetic/Data_in/LENGTH420CatExpe33SOPAS.Ptrim.trimmo.Maxee1.LABEL.sorted420NOPuntoComaBUE.MOD2.fas --db ../genetic/AlignedSeq/Hymenoptera_161_33SOPAS_UNIQUES_ExpeRelabel420_20_26_1_alignment.fasta --otutabout ../genetic/Data_out/ASVs_Table/allele_table_Hymenoptera_161_33SOPAS_UNIQUES_ExpeRelabel420_20_26_1_alignment.txt 

