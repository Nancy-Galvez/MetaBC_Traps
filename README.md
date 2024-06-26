**Review**
#
The repository contains the pipeline to perform bioinformatics tools from metabarcoding data of the project *Towards DNA metabarcoding-based haplotype for monitoring terrestrial arthropods communities*.
#
**Primers and overhang adapter:**

**COI** expected size 418 pb

B_F 5' CCIGAYATRGCITTYCCICG 3' (Shokralla et al., 2015)

Fol-degen-R 5’ TANACYTCNGGRTGNCCRAARAAYCA 3' (Yu et al., 2012)

The overhang adapter sequence must be added to the locus‐specific primer
for the region to be targeted. The Illumina overhang adapter sequences to be
added to locus‐specific sequences are:

Forward overhang: 5’ TCGTCGGCAGCGTCAGATGTGTATAAGAGACAG‐[locus‐
specific sequence]
Reverse overhang: 5’ GTCTCGTGGGCTCGGAGATGTGTATAAGAGACAG‐[locus‐
specific sequence]
#

## Repository organization
Contains data and scripts for the sections *Bioinformatics processing to identify OTUs at different thresholds of genetic similarity*, and *Community diversity and composition* of the manuscript.

The **bin** directory cointains the scripts used in this pipeline. All scripts run using `bin` as working directory. The **data** directory is not included in this repository, but data is available at this [dryad](https://XXXXXXX).

Data comes from Cornell Institute of Biotechnology, Cornell University, USA, for sequencing on a lane of Illumina MiSeq 2x300 bp.
#

### `/bin/`

The scripts in `/bin` should be run in the order they are numbered. R functions used by some of these scripts are not numbered and have the extension `.R`. Html notebooks are provided for some of the analyses in R.

Scripts content:

* `0.0install_software.txt` well, not actually a script, but this fille contains packages for the `0.1Processing_33BioSoup_mbc.sh` and `0.2Processing_33BioSoup_metaMATE.sh` scripts into Arribas et al., 2020. The packages are *fastqc*, *fastx-toolkit*, *trimmomatic-0.36*, *pairfq-0.17*, *usearch-9.2*, and *usearch-10*.
* `0.1Processing_33BioSoup_mbc.sh` Paired-end reads of samples were quality filtered following procedures described by Arribas et al. (2020). Briefly, processing included quality checking, primer removal, pair merging, quality filtering, denoising, and clustering each library independently. Size COI = 418 bp. 
* `0.2Processing_33BioSoup_metaMATE.sh` Paired-end reads of samples were quality filtered following procedures described by Arribas et al. (2020).  Briefly, processing included quality checking, primer removal, pair merging, quality filtering, denoising, and clustering each library independently. Size COI = 416-420 bp. 
* `0.3Steps_afterProcessing.txt`well, not actually a script. Contains the steps for each library: Get unix, blast to MEGAN and visualised tree in figtree.
* `0.4Searching_Stop_Codons.txt` not actually a script. Here, each ASV dataset was aligned in Geneious for searching codon stops. 
* `0.5to_get_ASV_tables.sh` script for geting a community table then generated with read-counts (haplotype abundance) of each retained ASV for the eight orders by matching ASVs against the complete collection of reads.
* `0.6to_get_UPGMAtree_GMYC_MH_lineages_All.r` This script gets all scripts for the first step (source) that obtained lineages at different clustering levels for each of the orders. STEP 1 and 2: We get analysis of the UPGMAtree and GMYC each lineages. Step 3: We get analysis to apply NODE.MIN get trees multiple each lineages (e.g. We went at "Arachnida/SpeciesDelimitation/2to apply NODE.MIN_Get_trees_multiple_each_lineage_ArachnidaStep2" and we ejecuted each script with differente order). 
* `2.MergertwoFiles_AllelesOtus&LimitesSpecies.txt` not actually a script. Here, each ASV dataset was joined with Species Limites in excel, after for searching "unoised". (AQUI).
* `3.to_get_conservative_threhold_on_original_all.r` This script gets conservative thershold on origibale table ASVs by artropods order. It is important to place the correct number of columns and rows, because this can change in each group. In this step, we can join all the groups or put them separately as in this case. (FALTA)
* `4.to_get_Diversity_All.r`This script gets all scripts (source) of diversity using 8 arthropods order at multi-hierarchical levels. (FALTA)
* `5.Plot_Diversity_All.r`This script plot "Community diversity and composition", "plot global Richness by sites", "Beta diversity", "Non-Metric Multidimensional scaling (NMDS) ordinations of community similarity", and "Accumulation Curves". (FALTA)
#

These scripts use the data in `genetic`, `figures` and `meta`.

### `/genetic`

Contains genetic *data in* and *data out* for each order

Genetic *data in* corresponds to `0.1Processing_33BioSoup_mbc.sh` output using the subset of each order. 


Genetic *data out* corresponds to `4.to_get_Diversity_All.r`using the subset of 5 orders community diversity and composition at the haplotype, 3% and 5% lineages for each of the five taxonomic orders studied.
#

### `/figures`

Contais figures of results:

* `/figures` from`5.Plot_Diversity_All.r` output using the subset of five orders for richness.

 `/figures` from`5.Plot_Diversity_All.r` output using the subset of Diptera and Hymenoptera order of the composition at the haplotype, 3% and 5% lineages.

### `/meta`
`ConservationForestNevadoToluca.csv` contains metadata for each of the samples sequenced in a lane Miseq. Each column names refer to:

* `id`: sample number ID 
* `code`: ID of the sequencing run 
 * `label_metabarcoding`: sample name of each library: e.g. `T_W_SJH_3SHCO5COW`: `T` Treatent), `W`Trap type, `SJH`San Juan de las Huertas (locality), `3SHCO5COW` ID of the sequencing sample run.
* `Short Description`: Description of the each soup.
* `Long Description`: More information for short description.
* `Number of read`: Total for each soup.
* `Merge_Pair`: Merge Pair for after filtering and denoising.
* `forest_type`: Name of the tree specie in the forest.

**END**
