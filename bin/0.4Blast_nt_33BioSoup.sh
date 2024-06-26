#!/bin/bash

#SBATCH -w keri
#SBATCH -n 10


#BLAST TO MEGAN

blastn -db /LUSTRE/Genetica/ngalvez/1SOPAS_Artropodos/4BLAST/blastdb_2021_12_20/NT_Nancy -query all.ZOTUS_UNIQUES_420_51ExpSINPuntoyComaFinal.fasta -outfmt 5 -out nt_USall.ZOTUS420_UNIQUES_51Conser.fas.xml -num_threads=4 -evalue 0.001 -max_target_seqs 100 &



