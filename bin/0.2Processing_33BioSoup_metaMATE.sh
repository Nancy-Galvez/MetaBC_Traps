#!/bin/bash

#SBATCH -w keri
#SBATCH -n 10

#1.Fastqc on all files

#fastqc *.fastq &
#(a) Descomprimir TODOS los archivos.

#$ for f in *.gz ; do gunzip -c "$f" > DATOS/"${f%.*}" ; done

#2.Descomprimir all files

#gzip -d *gz &

#3.PRIMER TRIMMING (fragment length 418)
cd ../genetic/raw
mkdir -p PTRIM416_420
for file in *R1.fastq
do
	fastx_trimmer -f 21 -i ${file} -o ./PTRIM416_420/${file%.*}.Ptrim.fq
done

for file2 in *R2.fastq
do
	fastx_trimmer -f 27 -i ${file2} -o ./PTRIM416_420/${file2%.*}.Ptrim.fq 
done


#4. FIRST QUALITY CONTROL
cd ./PTRIM416_420
mkdir -p TRIMMO

for file in *Ptrim.fq
do
       	java -jar ../../.././bin/TRIMMOMATIC/trimmomatic-0.36.jar SE -threads 20 -phred33 -trimlog ./TRIMMO/${file%.*}.log ${file} ./TRIMMO/${file%.*}.trimmo.fq TRAILING:20 
done


#PAIRING
cd ./TRIMMO
mkdir -p PAIRED 
for f1 in *R1.Ptrim.trimmo.fq
do
	echo $f1
	part1="$(echo $f1 | sed 's/\(.*\)R1.*/\1/')"
	f2="$part1"
	f2+="R2.Ptrim.trimmo.fq"
	echo "file1: $f1"
	echo "file2: $f2"
	echo "executing pairfq_lite ..."

	perl ../../../.././bin/pairfq_lite.pl makepairs -f $f1 -r $f2 -fp ./PAIRED/$f1.p.fq -rp ./PAIRED/$f2.p.fq -fs ./PAIRED/short2.1.s.fq -rs ./PAIRED/short2.2.s.fq
	echo "done"
done

#fastq_mergepairs
cd ./PAIRED
mkdir -p MergePair
for f1 in *R1.Ptrim.trimmo.fq.p.fq
do
	echo $f1
	part1="$(echo $f1 | sed 's/\(.*\)R1.*/\1/')" 

	f2="$part1"
        f2+="R2.Ptrim.trimmo.fq.p.fq"
	f3="$part1"
	
	f3+="R1_R2.Ptrim.trimmo.fastq"

	f4="$part1"
	
	f4+="R1_R2.Ptrim.trimmo.report.log"
	
	echo "file1: $f1"
	echo "file2: $f2"
	echo "file3: $f3"
	
	echo "executing fastq_mergepairs ..."
	
	../../../../.././bin/usearch9.2.64_i86linux32 -fastq_mergepairs $f1 -reverse $f2 -fastq_minovlen 50 -fastq_maxdiffs 15 -fastq_maxdiffpct 50 -fastqout ./MergePair/$f3 -report ./MergePair/$f4 
	echo "done"
done

#MAXEE1

cd ./MergePair
mkdir -p maxee1
for f1 in *R1_R2.Ptrim.trimmo.fastq
do
	echo $f1
	part1="$(echo $f1 | sed 's/\(.*\)fastq*/\1/')"    
	f2="$part1"
        f2+="Maxee1.fas"

	f3="$part1"
        f3+="Maxee1.log"
	
	echo "input: $f1"
	echo "output: $f2"
	echo "logfile: $f3"

	echo "executing maxee 1 ..."
	
	../../../../../.././bin/usearch9.2.64_i86linux32 -fastq_filter $f1 -fastq_maxee 1 -fastaout ./maxee1/$f2 &> ./maxee1/$f3
	echo "done"
done

#LABELLING
cd ./maxee1
mkdir -p LABEL
for file in *.fas
do
	
	part1="$(echo $file | sed 's/\(.*\)_R1.*.*/\1/')" 
	#part2="$(echo $file | sed 's/.*_001\(.*\).Ptrim.*/\1/')"
	
	label="$part1"		
	echo $file
	echo $label
	sed "-es/^>\(.*\)/>\1;barcodelabel=$label;/" < ${file} > ./LABEL/${file%.*}.LABEL.fas 
done


#sortbylength and keeping only reads of 416-420 (Amplicon expected length of 418, +- 2bp)

cd ./LABEL
mkdir -p LENGTH416-420
for f1 in *LABEL.fas
do	
	echo $f1
	part1="$(echo $f1 | sed 's/\(.*\)fas*/\1/')" 
	f2="$part1"
        f2+="sorted420.fas"
	f2log="$part1"
	f2log+="sorted420.log"
				
	echo "input: $f1"
	echo "output: $f2"
	echo "executing sortbylength 420 ..."
	
	../../../../../../../.././bin/usearch10.0.240_i86linux32 -sortbylength $f1  -fastaout ./LENGTH416-420/$f2 -minseqlength 416 -maxseqlength 420 &> ./LENGTH416-420/$f2log
	echo "done sortbylength 420 ..."
done

#get uniques by lib
cd ./LENGTH416-420
mkdir -p UNIQUESbyLIB
for f1 in *sorted420.fas
do	
	echo $f1
	part1="$(echo $f1 | sed 's/\(.*\)fas*/\1/')" 

	f2="$part1"
        f2+="uniques.fas"
	f2log="$part1"
	f2log+="uniques.log"

	echo "input: $f1"
	echo "output: $f2"			
	echo "executing fastx_uniques ..."
	../../../../../../../../.././bin/usearch10.0.240_i86linux32 -fastx_uniques ./$f1 -fastaout ./UNIQUESbyLIB/$f2 -sizeout &> ./UNIQUESbyLIB/$f2log
	echo "done fastx_uniques ..."	
done				
				
#unosie3 by lib
cd ./UNIQUESbyLIB
mkdir -p UNOISE
for f1 in *uniques.fas
do
	echo $f1
	part1="$(echo $f1 | sed 's/\(.*\)fas*/\1/')" 			
	f2="$part1"
        f2+="zotus.fas"
	f2log="$part1"
	f2log+="zotus.log"
	f2tabout="$part1"
	f2tabout+="zotus.tabbedout.txt"
						
	echo "input: $f1"
	echo "output: $f2"
	echo "executing unoise3 ..."
	../../../../../../../../../.././bin/usearch10.0.240_i86linux32 -unoise3 ./$f1 -zotus ./UNOISE/$f2 -minsize 4 -tabbedout ./UNOISE/$f2tabout &>./UNOISE/$f2log
	echo "done unoise3 ..."
done

#labelling_zotus.sh
 
cd ./UNOISE
mkdir -p ./RELABEL
for file in *zotus.fas
do
	part1="$(echo $file | sed 's/\(.*\)_R1.*.*/\1/')" 
	#part2="$(echo $file | sed 's/.*_001\(.*\).Ptrim.*/\1/')"
	
	label="$part1"	
	
	echo $file
	echo $label
	sed "-es/^>\(.*\)/>\1_$label/" < ${file} > ./RELABEL/${file%.*}.relabel.fas
done

#juntar todos los zotus individuales

#cd ./RELABEL
#cat *relabel.fas > all.ZOTUSbySAMPLE.ZONA.fas


#Dejar solo los únicos
#usearch100 -fastx_uniques all.ZOTUSbySAMPLE.ZONA.fas -fastaout all.ZOTUS_UNIQUES.ZONA.fas -sizeout &> all.ZOTUS_UNIQUES.ZONA.log

#BLAST TO MEGAN
#blastn -db /db/Paula/nt_Sanger -query all.ZOTUS_UNIQUES.ZONA.fas -outfmt 5 -out nt_SangerVsall.ZOTUS_UNIQUES.ZONA.fas.xml -num_threads=4 -evalue 0.001 -max_target_seqs 100 &



#############NOTAS SOBRE METAMATE######################################################################################################
#NOTA1: tras hacer el "Blast para MEGAN", se usa MEGAN para visualizar y exportar los set de datos de cada uno de los 8 órdenes. Los set de datos de ahora serán más grandes,
#ya que incluyen secuencias con +- 2 bp sobre la longitud esperada del amplicón. Estos 8 set de datos, se renombran para que el nombre de cada ASV incluya el orden taxonómico
# al que pertenece (en Geneious o con comando sed), y se concatenan para servir como imput de ASV para metaMATE.

#NOTA2:El otro input para metaMATE será el archivo con todos los reads (nombrados por librería, ver indicaciones en Tutorial metaMATE) antes de hacer el "get uniques", es decir,
#el que se usa para asignar la abundancia a read-count con que aparece cada ASV. En este caso, hay que usar los nuevos archivos generados para cada librería con 416-420 bp,
#renombrados por librería y concatenados como un solo archivo (comprobad esto último en el tutorial de metaMATE)

#NOTA3:El tercer y último input para metaMATE son las sequencias de refrencia, en este caso todos los bc de BOLD para cada Orden.
###################################################################################################################



#juntar todos los reads de 420
#cd /db/Paula/PRUEBAS/P2/PTRIM/TRIMMO/PAIRED/MergePair/maxee1/LABEL/LENGTH420
#cat *sorted420.fas > /db/Paula/PRUEBAS/P2/PTRIM/TRIMMO/PAIRED/MergePair/maxee1/LABEL/LENGTH420/UNIQUESbyLIB/UNOISE/RELABEL/ALLreads.ZONA.420.fas


#final step after preparing file in geneious
#usearch100 -search_exact all.420.fas  -db All_s2_sinSTOPS_revised.fasta -strand plus -otutabout allele_table_All_s2_sinSTOPS_revised.txt  
