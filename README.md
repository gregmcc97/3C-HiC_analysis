# 3C-HiC_analysis manuscript

Code used for manuscript methods

## Software used:
- [PrinSeq-lite](http://prinseq.sourceforge.net/manual.html#STANDALONE) v0.20.4
- [Cutadapt](https://github.com/marcelm/cutadapt) v2.5
- [Bowtie2](https://github.com/BenLangmead/bowtie2) v2.3.4.1
- [SAMtools](https://github.com/samtools/samtools) v0.1.19
- [BEDTools](https://github.com/arq5x/bedtools2) v2.25.0
- [MetaPhlAn](https://github.com/biobakery/MetaPhlAn) v3.0
- [Megahit](https://github.com/voutcn/megahit) v1.1.3
- [ABRicate](https://github.com/tseemann/abricate) v0.9.8
- [CoverM](https://github.com/wwood/CoverM) v0.4.0
- [Burrow-Wheeler Aligner](https://github.com/lh3/bwa) v0.7.12
- [Prokka](https://github.com/tseemann/prokka) v1.14.6
- [Kraken2](https://github.com/DerrickWood/kraken2) v2.0.8
  - using [kraken2-microbial database](https://lomanlab.github.io/mockcommunity/mc_databases.html)
- [BLAST command line](https://www.ncbi.nlm.nih.gov/books/NBK279690/) v2.2.31
- [GNU Parallel](https://www.gnu.org/software/parallel/sphinx.html) v20190322

For making heatmaps:
- [R](https://www.r-project.org/) v4.0.0 - with packages:
  - [pheatmap](https://github.com/raivokolde/pheatmap)
  - [RColorBrewer](https://cran.r-project.org/web/packages/RColorBrewer/index.html)

## Processing reads
Reads were first deduplicated using PrinSeq-lite:
```
for i in *_1.fastq ; do prinseq-lite.pl -fastq $i -fastq2 ${i/_1.fastq/_2.fastq} -derep 14 ; done
for i in *_1_prinseq_good* ; do mv $i ${i%_1_prinseq_good_*}_derep_1.fastq ; done
for i in *_2_prinseq_good* ; do mv $i ${i%_2_prinseq_good_*}_derep_2.fastq ; done
```
Adapters were trimmed and low-quality reads were filtered using Cutadapt with a fasta file of adapters [ADAPTER_FILE].fasta. Note that `--nextseq-trim=20` was used in below command - this option was used if library was sequenced on an instrument that uses two-color chemistry (Illumina NextSeq or NovaSeq, see [CutAdapt manual](https://cutadapt.readthedocs.io/en/stable/guide.html) for more information). If not, then `-q 20` was used as the quality filtering option instead.
```
for i in *derep_1.fastq ; do cutadapt -j 16 --nextseq-trim=20 -b file:[ADAPTER_FILE].fasta -B file:[ADAPTER_FILE].fasta -m 60 -o ${i/derep_1.fastq/trimmed_1.fastq} -p ${i/derep_1.fastq/trimmed_2.fastq} $i ${i/_1.fastq/_2.fastq} ; done
```
Host DNA was removed following [this tutorial](https://www.metagenomics.wiki/tools/short-read/remove-host-sequences).

For human faecal samples, [GRCh38.p13 human reference genome](https://www.ncbi.nlm.nih.gov/assembly/GCF_000001405.39/) was downloaded and a Bowtie2 index called human_DB was created. For mouse faecal samples, [GRCm38.p6 mouse reference genome](https://www.ncbi.nlm.nih.gov/assembly/GCF_000001635.26/) was used. Then:
```
for i in *trimmed_1.fastq ; do bowtie2 --threads 16 -x human_DB -1 $i -2 ${i/_1.fastq/_2.fastq} -S ${i/_trimmed_1.fastq/_mapped_and_unmapped.sam} ; done 
for i in *_mapped_and_unmapped.sam ; do samtools view -@ 16 -bS $i > ${i/.sam/.bam} ; done
for i in *_mapped_and_unmapped.bam ; do samtools view -@ 16 -b -f 12 -F 256 $i > ${i/_mapped_and_unmapped.bam/_bothEndsUnmapped.bam} ; done
for i in *_bothEndsUnmapped.bam ; do samtools sort -@ 16 -n $i > ${i/.bam/_sorted.bam} ; done
for i in *_bothEndsUnmapped_sorted.bam ; do bedtools bamtofastq -i $i -fq ${i/_bothEndsUnmapped_sorted.bam/_hr_1.fastq} -fq2 ${i/_bothEndsUnmapped_sorted.bam/_hr_2.fastq} ; done
```

## Taxonomic profiling of reads
Reads were taxonomically profiled using MetaPhlAn v3.0:
```
mkdir metaphlan_output
#including unknown estimation:
for i in *_hr_1.fastq ; do metaphlan ${i},${i/_1.fastq/_2.fastq} --bowtie2out metaphlan_output/metaphlan_${i/_hr_1.fastq/}.bowtie2.bz2 --nproc 16 --input_type fastq --unknown_estimation --add_viruses > metaphlan_output/${i/_hr_1.fastq/_profiled_metagenome.txt} ; done
#classification of known reads only:
cd metaphlan
for i in metaphlan*.bz2 ; do metaphlan $i --nproc 16 --input_type bowtie2out --add_viruses -o ${i/.bowtie2.bz2/_knownprofile.txt} ; done
```
A table for creating a figure of the class abundance was generated following [this python script](https://github.com/flannsmith/metaphlan-plot-by-taxa/blob/72c0db5302e8b1e05732ef46a1252cad72ff8075/Converting%20Metaphlan%20profile%20to%20Phyloseq%20objects.ipynb) after a merged table of all results was created using MetaPhlAn's `merge_metaplan_tables.py` script:
```
merge_metaphlan_tables.py *_profiled_metagenome_known.txt > merged_abundance_table_known.txt

#table adjusted to remove second column
cut -f 1,3-200 merged_abundance_table_known.txt > merged_abundance_table_known2.txt
#column titled clade_name was then manually renamed to ID before python script was used
```
[GraphPad Prism](https://www.graphpad.com/scientific-software/prism/) was used to create the figure.

## Metagenomic assembly and identification of antimicrobial resistance gene contigs
Reads were assembled using Megahit. For Hi-C datasets, accompanying shotgun metagenomic reads were used for assembly. For meta3C datasets, the 3C reads were used for assembly.
```
for i in *_hr_1.fastq ; do megahit -t 16 --min-contig-len 1000 --k-min 21 --k-max 141 --k-step 12 --merge-level 20,0.95 -1 $i -2 ${i/_1.fastq/_2.fastq} -o megahit_assembled_${i/_hr_1.fastq/} ; done
```
Contigs containing antimicrobial resistance genes (ARGs) were then identified:
```
for i in *_final.contigs.fa ; do abricate --db resfinder --threads 16 --minid 95 --mincov 75 $i > ${i/.fa/.abricate.tsv} ; done
```
The relative abundance of the ARGs was then calculated.

First, the ARG sequences were isolated from the ARG contigs:
```
#in directory containing assemblies directories
for i in megahit_assembled_* ; do mkdir ${i}/${i/megahit_assembled_/}_arg_contigs ; done

for i in megahit_assembled_* ; do grep -v "SEQUENCE" ${i}/${i/megahit_assembled_}_final.contigs.abricate.tsv | while read sample contig start end strand arg_full rest ; do arg=$(echo "$arg_ful" | cut -f6 | sed s/"'"/""/g | sed s/"-"/""/g | sed s/"\."/""/g | sed 's/[()]//g') && size=$(echo "${end} - ${start} + 1" | bc) && length=$(grep "\<${contig}\>" ${i}/${i/megahit_assembled_}_final.contigs.fa | cut -d" " -f4 | sed 's!len=!!') && distance=$(echo "${length} - ${end} + 1" | bc) && sequence=$(grep -A1 "\<${contig}\>" ${i}/${i/megahit_assembled_}_final.contigs.fa | awk 'BEGIN{ cmd = "cut -c '${start}'- | rev | cut -c '${distance}'- | rev" } NR % 2 { print } NR % 2 == 0 { print | cmd; close(cmd) }' | tail -1) && title=$(grep "\<${contig}\>" ${i}/${i/megahit_assembled_}_final.contigs.fa | cut -d" " -f1) && (echo ${title}_${arg}_len=${size} && echo $sequence) >> ${i}/${i/megahit_assembled_}_arg_contigs/${i/megahit_assembled_}_all_args.fa ; done ; done
```
Then run CoverM was used to get counts of reads (3C/Hi-C reads and respective shotgun reads) that align to each ARG sequence for each sample:
```
coverm contig -r [SAMPLE]_all_args.fa -t 16 --coupled [read_1].fastq [read_2].fastq -m count > [SAMPLE]_all_args_counts.txt
```
Then use counts to calulcate RPKM (where [READ_COUNTS_TABLE].tsv is a table of total number of reads for each set of reads):
```
for i in megahit_assembled_* ; do echo -e CONTIG'\t'START'\t'END'\t'ARG_name'\t'LENGTH'\t'READS_MAPPED'\t'TOTAL_READS'\t'RPKM'\t'LOG_RPKM > ${i/megahit_assembled_/}_arg_rpkm.tsv && grep -v "SEQUENCE" ${i}/${i/megahit_assembled_}_final.contigs.abricate.tsv | while read sample contig start end strand arg_full rest ; do arg=$(echo "$arg_ful" | cut -f6 | sed s/"'"/""/g | sed s/"-"/""/g | sed s/"\."/""/g | sed 's/[()]//g') && size=$(echo "${end} - ${start} + 1" | bc) && length=$(grep "\<${contig}\>" ${i}/${i/megahit_assembled_}_final.contigs.fa | cut -d" " -f4 | sed 's!len=!!') && count=$(grep "${contig}_${arg}" [SAMPLE]_all_args_counts.txt | cut -f2) && && total=$(grep "${i/megahit_assembled_/}_[TYPE_OF_READ]" [READ_COUNTS_TABLE].tsv | cut -f2 | tr -d $'\r') && rpkm=$(echo "scale=10 ; ( ${count} / ( ${total} / 1000000 ) ) / ( ${length} / 1000 )" | bc | awk '{printf "%.10f\n", $0}') && logrpkm=$(echo "l(${rpkm})/l(10)" | bc -l | awk '{printf "%.10f\n", $0}') && echo -e ${contig}'\t'${start}'\t'${end}'\t'${gene}'\t'${length}'\t'${count}'\t'${total}'\t'${rpkm}'\t'${logrpkm} >> ${i/megahit_assembled_/}_arg_rpkm.tsv ; done ; done
```


IS elements were also identified using [ISfinder database](https://github.com/thanhleviet/ISfinder-sequences) (make custom ABRicate database)
```
cd [directory/containing/abricate/db]
mkdir is_seq
wget https://github.com/thanhleviet/ISfinder-sequences/raw/master/IS.fna
mv IS.fna sequences
abricate --setupdb

#then return to assembly to run abricate:
for i in *.fa ; do abricate --db is_seq --threads 16 --minid 99 --mincov 60 $i > ${i/.fa/[ABRICATE_IS].tsv} ; done
```

## Mapping 3C/Hi-C reads and finding intercontig reads
3C/Hi-C reads were mapped to their respective metagenomic assembly with [bwa](https://github.com/lh3/bwa) using the aln and sampe commands.

To map only the first 50 bp of the Hi-C reads, the reads were first trimmed to keep only the first 50 bp:
```
for i in *.fastq ; do cut -c 1-50 $i > ${i/_hr_/_hr_cut_} ; done
```
Bwa indexes were generated:
```
mkdir bwa_index
cd bwa_index
bwa index -p [index_name] -a is [directory/containing/assembly]/[ASSEMBLY].fa
```
Reads were mapped & filtered for reads that map with mapping quality > 20:
```
mkdir mapped
for i in *_hr_cut_*.fastq ; do bwa aln -t 16 bwa_index/[index_name] $i > mapped/${i/.fastq/.aligned.sai} ; done

#create sam file from alignment files
for i in *_hr_cut_1.fastq ; do bwa sampe bwa_index/[index_name] mapped/${i/.fastq/.aligned.sai} mapped/${i/_1.fastq/_2.aligned.sai} $i ${i/_1.fastq/_2.fastq} > mapped/${i/_1.fastq/_aligned.sam} ; done
cd mapped

#convert to bam file (can now delete sam file if space needed)
for i in *aligned.sam ; do samtools view -@ 16 -bSh $i > ${i/.sam/.bam} ; done

#filter alignment file to get reads that map with mapping quality >20
for i in *aligned.bam ; do samtools view -@ 16 -bh -q 20 > ${i/.bam/_filtered.bam} ; done

#if you want stats of how many reads mapped:
for i in *.bam ; do samtools flagstat -@16 $i > ${i/.bam/_flagstat.txt} ; done
```
Using the filtered alignment, intercontig reads (3C/Hi-C reads where each read of the pair maps to a different contig) were isolated:
```
for i in *_aligned_filtered.bam ; do samtools view -F 14 $i | awk '$7!="=" {print $0}' > ${i/_aligned_filtered.bam/_intercontig.sam} ; done

#samtools view -F 14 filters out reads that are unmapped, have an unmapped mate, or are mapped in a proper pair
# grep -v "=" removes remaining aligned reads that align to the same contig as their mate ('=' in RNEXT column of sam file (column 7))
```
The same mapping process was carried out using the shotgun metagenomic reads (if the sample had accompanying shotgun reads) to find spurious intercontig reads in the shotgun datasets.

## Analysis of G_3C reads mapping to spike-ins 
Complete genome sequences of the [E. coli E3090](https://www.ncbi.nlm.nih.gov/biosample/SAMEA4699317) and [E. faecium E745](https://www.ncbi.nlm.nih.gov/biosample/SAMN04045274/) spike-ins were downloaded from NCBI. Whole-genome sequence reads and 3C reads were mapped to the spike-in genomes as described above. 

To see what regions of the spike-in genomes the reads were mapping to, the genomes were annotated using [Prokka](https://github.com/tseemann/prokka):
```
prokka --outdir prokka_output_[SPIKE-IN] --prefix [SPIKE-IN] --cpus 4 --addgenes [SPIKE-IN].fasta 
```
The relevant information from the .gff output files was then added to a table:
```
# [Genome] = corresponding spike-in genome
grep "CDS" [genome].gff| sed -e 's/;/\t/g' |cut -f 4,5,9,13 > [Genome]_genes.tsv
# grep "CDS" e745.gff # Strips lines containing "CDS", ignoring headers and sequence.
# | sed -e 's/;/\t/g' # Replaces ; with a tab, to allow separation of important info 
# |cut -f 4,5,9,13 # Uses tabs to differentiate between fields, prints start, stop, Locus tag and product.
```
The sam header was also isolated:
```
tail -n +9 [mapping file of either WGS or 3C reads].sam > [Genome]_sam_lookup.tsv
```
These files were then used in R to find the mapped genome regions:
```
# In R
# Install and load packages:
install.packages("dplyr")
install.packages("tidyr")
install.packages("purrr")
library(dplyr)
library(tidyr)
library(purrr)
# Set working directory
setwd('/[working_directory]')
# Read file generated from .gff
[Genome]_genes <- read.csv(file="[Genome]_genes.tsv",header=FALSE,sep='\t')
# Set column names
colnames([Genome]_genes) <- c('Start','Stop','Locus', 'Product')
# Create new df with seq 1:max(GENOME_LENGTH)
# Can check max coord position with max([Genome]_genes$Stop)
lookup_[Genome] <- data.frame(Pos = seq(1,[GENOME_LENGTH]))
# Match columns based on start/stop co-ords
matched_[Genome] <- [Genome]_genes %>% mutate(Pos = map2(Start, Stop, `:`)) %>%
  unnest(Pos) %>% select(3:5) %>% right_join(lookup_[Genome]) %>%
  arrange(Pos) %>% select(3,1,2)
# Read processed sam file
[Genome]_sam <- read.csv(file='[Genome]_sam_lookup.tsv', header=FALSE, sep='\t')
# Extract interesting columns
sam_lookup_[Genome] <- [Genome]_sam[,c(1,4)]
# Set column names
colnames(sam_lookup_[Genome]) <- c('Info', 'Pos')
# Match mapped reads to genes
output_[Genome] <- merge(x = matched_[Genome], y = sam_lookup_[Genome], by = 'Pos')
write.table(output_[Genome], file = "[Genome]_[WGS/3C]_links.tsv", quote=FALSE, sep='\t', row.names = FALSE)
```
The output table was then used to assign mapped genome regions using the following bash script:
```
for i in *links.tsv ; do cat $i | while read line ; do product=$(echo "$line" | cut -f3) && if [[ $product == "NA" ]] ; then echo -e "${line}"'\t'intergenic ; else if [[ $product == *"IS"* ]] ; then echo -e "${line}"'\t'is_element ; else if [[ $product == *"transposase"* ]] ; then echo -e "${line}"'\t'transposon ; else if [[ $product == *"hypothetical protein"* ]] ; then echo -e "${line}"'\t'hypotehtical_protein ; else if [[ $product == *"prediction"* ]] ; then echo -e "${line}"'\t'hypotehtical_protein ; else if [[ $product == *"product=putative protein"* ]] ; then echo -e "${line}"'\t'hypotehtical_protein ; else if [[ $product == *"gene="* ]] ; then echo -e "${line}"'\t'annotated_gene ; else if [[ $product == *"locus_tag"* ]] ; then echo -e "${line}"'\t'annotated_gene ; else if [[ $product == *"db_xref"* ]] ; then echo -e "${line}"'\t'annotated_gene ; else if [[ $product == *"protein"* ]] ; then echo -e "${line}"'\t'annotated_gene ; else if [[ $product == *"note"* ]] ; then echo -e "${line}"'\t'annotated_gene ; else if [[ $product == *"product"* ]] ; then echo -e "${line}"'\t'annotated_gene ; else echo -e "${line}"'\t'other ; fi ; fi ; fi ; fi ; fi ; fi ; fi ; fi ; fi ; fi ; fi ; fi ; done > ${i/.tsv/_labelled.tsv} ; done
```
## Filtering intercontig reads
To minimise problematic noise from spurious intercontig reads (intercontig reads not originating from cross-linked fragments of DNA), intercontig reads that map within the first or last 500 nt of a contig are filtered out:
```
# First find lengths of contigs in assembly:
grep ">" [ASSEMBLY].fa | cut -d" " -f1,4 | sed 's!len=!!' | sed 's/\s/\t/'g > [ASSEMBLY]_contig_lengths.tsv

# Output seperate sam files containing reads mapping <500 nt (within) and >500 nt (not) of the start of end of a contig:
for i in *intercontig.sam ; do cat $i | while read -r name qual mappedto position rest ; do grep "\<${mappedto}\>" [ASSEMBLY]_contig_lengths.tsv | while read -r contig length ; do if (( $position < 501 )) ; then echo -e ${name}'\t'${qual}'\t'${mappedto}'\t'${length}'\t'${position}'\t'"${rest}" >> ${i/.sam/_within.sam} ; else if distance=$(echo "${length} - ${position}" | bc) && (( $distance < 501 )) ; then echo -e ${name}'\t'${qual}'\t'${mappedto}'\t'${length}'\t'${position}'\t'"${rest}" >> ${i/.sam/_within.sam} ; else echo -e ${name}'\t'${qual}'\t'${mappedto}'\t'${length}'\t'${position}'\t'"${rest}" >> ${i/.sam/_not.sam} ; fi ; fi ; done ; done ; done
```
To calculate the proportion of intercontig reads mapping <500 nt from ends of contig:
```
# Prepare cut first, third, and fourth column from sam file
for i in *intercontig.sam ; do cut -f1,3,4 $i > ${i/.sam/.tsv} ; done
# Calculate position of alignment for all reads
for i in *intercontig.tsv ; do cat $i | while read -r name mappedto position ; do grep "\<$mappedto\>" [Assembly]_contig_lengths.tsv | while read -r contig length ; do if (( $position < 501 )) ; then echo -e ${name}'\t'${mappedto}'\t'${length}'\t'${position}'\t'${position}'\t'within'\t'start ; else if distance=$(echo "${length} - ${position}" | bc) && (( $distance < 501 )) ; then echo -e ${name}'\t'${mappedto}'\t'${length}'\t'${position}'\t'${distance}'\t'within'\t'end ; else echo -e ${name}'\t'${mappedto}'\t'${length}'\t'${position}'\t'${distance}'\t'not'\t' ; fi ; fi ; done ; done >  ${i/.tsv/_positions.tsv} ; done
# Find total number of reads mapping within 500 nt of the ends of a contig or not 
# File types.txt = file containing lines "within" and "not"
for i in *positions.tsv ; do cat types.txt| while read line ; do grep "\<$line\>" $i | wc -l | cat | while read word ; do echo -e ${word}'\t'${line} ; done ; done > ${i/.tsv/_totals.tsv} ; done
```
This was also performed on non-intercontig reads to find the proportion mapping <500 nt from ends of contig. Non-intercontig reads were isolated using:
```
#remove unmapped reads, reads with unmapped mate. Then keep only reads mapped in a proper pair:
for i in *aligned_filtered.bam ; do echo "samtools view -F12 -f2 -@4 $i > ${i/.bam/_nonintercontig.sam}" ; done
```
## Link ARGs to hosts 
ARGs were then linked to thier hosts using the 3C/Hi-C intercontig reads. Several data files were first set up in a single directory (calling it [DATA/DIRECTORY] here) for ease of access:
- [SAMPLE]_arg_contigs.tsv - two columns showing contig_name,ARG_name (if multiple ARGs were on one contig, the contig is listed once with the ARG names in column 2 merged e.g. CONTIG_NAME  mph(E)_1-msr(E)_1):

`cut -f2,6 [ABRICATE_ARG].tsv | grep -v "SEQUENCE" > [DATA/DIRECTORY]/[SAMPLE]_arg_contigs.tsv`
- [SAMPLE]_is_list - file containing IS element contig names only:

`cut -f2 [ABRICATE_IS].tsv | grep -v "SEQUENCE" > [DATA/DIRECTORY]/[SAMPLE]_is_list`
- [SAMPLE]_contigs.fa - assembly file:

`cp [ASSEMBLY.fa] [DATA/DIRECTORY]/[SAMPLE]_contigs.fa`


The following script was then used to link ARGs to their hosts. For classifying the hosts, BLASTN is used to identify plasmid contigs (using [NCBI's nucleotide (nt) database](https://ftp.ncbi.nlm.nih.gov/blast/db/)). Kraken2 is then used to identify the hosts. For the Kraken2 database, the prebuilt [kraken2-microbial-database](https://lomanlab.github.io/mockcommunity/mc_databases.html) was downloaded and used.

```
#create variables to easily access all the data files:
export READ_SUFFIX=[SAMPLE]
export DATASET_DIR=[DATA/DIRECTORY]
export HEATMAP_DIR=[DIRECTORY/WHERE/YOU/WANT/YOUR/HEATMAP/TABLE/FILES]
export ARG_CONTIGS=${DATASET_DIR}/${READ_SUFFIX}_arg_contigs.tsv
export IS_LIST=${DATASET_DIR}/${READ_SUFFIX}_is_list
export ASSEMBLY=${DATASET_DIR}/${READ_SUFFIX}_contigs.fa

#ONLY APPLICABLE FOR SAMPLES WITH MULTIPLE EZNYMES USED FOR 3C/Hi-C (with separate reads):
export ENZYME_1=[enzyme 1 suffix e.g. _H]
export ENZYME_2=[enzyme 2 suffix e.g. _M]

#in mapped directory:
#find contigs linked to ARGs (individual file for each ARG, named as the ARG contig name e.g. k141_1234)
mkdir ${READ_SUFFIX}_linked_contigs 
mkdir ${READ_SUFFIX}_linked_contigs/cut
tab=$'\t'
for i in ${READ_SUFFIX}*_intercontig_not.sam ; do cut -f1 $ARG_CONTIGS | while read contig ; do grep "$contig$tab" $i > ${READ_SUFFIX}_linked_contigs/${contig}_${i/_intercontig_not.sam/} && cat ${READ_SUFFIX}_linked_contigs/${contig}_${i/_intercontig_not.sam/} | while read line ; do column=$(echo "$line" | awk -v b="${contig}" '{for (i=1;i<=NF;i++) { if ($i == b) { print i } }}') && if [[ "${column}" == "3" ]] ; then echo "$line" | cut -f1,9 >> ${READ_SUFFIX}_linked_contigs/cut/${contig}_${i/_intercontig_not.sam/} ; else if [[ "${column}" == "9" ]] ; then echo "$line" | cut -f1,3 >> ${READ_SUFFIX}_linked_contigs/cut/${contig}_${i/_intercontig_not.sam/} ; fi ; fi ; done ; done ; done
cd ${READ_SUFFIX}_linked_contigs/cut

#remove duplicates (currently have duplicate lines for both read 1 and 2 from the same pair)
mkdir unique
for i in k141_* ; do cat $i | sort | uniq > unique/$i ; done 
cd unique 

#get list of unique contigs linked to ARG 
mkdir uniq 
for i in k141_* ; do cut -f2 $i | sort | uniq > uniq/$i ; done 
cd uniq 

#get count for how many times each unique contig is linked to ARG
mkdir counts 
for i in k141_* ; do cat $i | while read line ; do grep "\<$line\>" ../$i | wc -l | cat | while read word ; do echo -e ${word}'\t'${line} >> counts/$i ; done ; done ; done 
cd counts 

#filter links so that a contig is only considered linked if linked by >=5 intercontig read pairs
#also removes linked IS element contigs 
mkdir unique_filtered 
for i in k141* ; do sort $i | uniq | sort -nr | grep -v -P '^1\tk141'\|'^2\tk141'\|'^3\tk141'\|'^4\tk141' | grep -vf $IS_LIST > unique_filtered/${i} ; done 
cd unique_filtered 

#get sequences for linked contigs
mkdir contigs
for i in k141_* ; do cat $i | while read count contig ; do grep -A1 "\<${contig}\>" $ASSEMBLY >> contigs/${i}_contigs.fa ; done ; done 

#classify linked contigs
cd contigs 
#first blast them:
mkdir blasted 
for i in k141* ; do blastn -query $i -db /DB/BLAST_DB/nt2/nt -num_threads 16 -max_target_seqs 5 -max_hsps 1 -outfmt "7 qseqid sacc bitscore evalue qcovs pident stitle" > blasted/${i/_contigs.fa/_blasted} ; done
#then classify
mkdir kraken 
for i in k141* ; do kraken2 --threads 16 --use-names --db /DB/KRAKEN_DB/kraken2-microbial $i --output kraken/${i/_contigs.fa/_kraken} ; done 
cd .. 
#get top BLAST results for each contig:
mkdir contigs/blast_first_results 
for i in k141* ; do cut -f2 $i | while read line ; do grep "\<$line\>" contigs/blasted/${i}_blasted | head -2 | tail -1 >> contigs/blast_first_results/${i} ; done ; done
mkdir pasted
for i in k141* ; do cat $i | while read count contig ; do blast=$(grep "\<$contig\>" contigs/blast_first_results/$i | cut -f7) && kraken=$(grep "\<$contig\>" contigs/kraken/${i}_kraken | cut -f3) && echo -e ${count}'\t'${contig}'\t'${blast}'\t'${kraken} >> pasted/$i ; done ; done
cd pasted

#start making heatmap file 
#get list of linked contig classifications and counts. At this point, contigs that BLAST identified as plasmid DNA are labelled as "Plasmid DNA", and the rest are labelled with their Kraken2 classification (down to genus level)
#also this is written for zsh - if using BASH, change ":u" to "^^" e.g. [[ "${classification^^}" == *"PLASMID"* ]]
mkdir names
for i in k141* ; 
do cat $i | while IFS=$'\t' read count contig blast kraken ; 
do if [[ "${blast:u}" == *"PLASMID"* ]] ; then
echo -e ${count}'\t'"Plasmid DNA" ;
else kraken_classification=$(echo "${kraken}" | cut -d" " -f1 | sed 's/(taxid//' | sed 's/[][]//g' | sed s/"'"/""/g | sed s/"-"/""/g | sed s/"\."/""/g) 
&& echo -e ${count}'\t'${kraken_classification} ;
fi >> names/$i ; done ; done 
cd names

#if multiple enzymes are used, now combine the files. If not then skip:
mkdir combined 
for i in k141*$ENZYME_1* ; do cat $i >> combined/${i/$ENZYME_1/} ; done
for i in k141*$ENZYME_2* ; do cat $i >> combined/${i/$ENZYME_2/} ; done
cd combined

#get proportions for links to each unique classification (contigs labelled "Plasmid DNA" removed at this point):
mkdir added 
for i in k141_* ; do total=$(grep -v "Plasmid DNA" $i | cut -f1 | paste -sd+ | bc) && grep -v "Plasmid DNA" $i | while read count name ; do combined=$(grep "\<${name}\>$" $i | cut -f1 | paste -sd+ | bc) && proportion=$(echo "scale=6 ; ${combined} / ${total}" | bc | awk '{printf "%.6f\n", $0}') && echo -e ${combined}'\t'${name}'\t'${proportion} ; done | sort | uniq > added/$i ; done
cd added

#get list of all classifications linked to ARGs
for i in k141_* ; do cat $i | cut -f2 ; done | sort | uniq > classification_list

#make heatmap table
#first convert each file into a list of all classifications and proportion of links to each (no links = 0, all links = 1)
mkdir columns 
for i in k141_* ; do cat classification_list | while read name ; do cat $i | while IFS=$'\t' read count title proportion ; do if [[ "${title}" == "${name}" ]] ; then echo -e ${name}'\t'${proportion} ; else echo -e ${name}'\t'0 ; fi ; done | sort -r | uniq | grep -m1 "\<${name}\>" ; done > columns/columns_${i} ; done
cd columns 

#add ARG names to top of list (the arg_rem part is so the filename contains the ARG without any special characters)
for i in columns_k141_* ; do b=${i/columns_/} && c=${b/_${READ_SUFFIX}/} && grep "\<${c}\>" $ARG_CONTIGS | while read contig arg ; do arg_rem=$(echo ${arg//[\(\)]/} | sed s/"'"/""/g | sed s/"-"/""/g | sed s/"\."/""/g) && (echo ${arg} && cut -f2 $i) > ${i/columns_/}_${arg_rem} ; done ; done

#get list of classifications linked to ARGs with blank first line (for column 1 of heatmap)
(echo -en '\n' && cat ../classification_list) > heatmap_list 

#make and run command for creating heatmap table, and copy the table to your heatmap output directory
(echo paste heatmap_list && for i in k141_* ; do echo $i ; done && echo "| sed 's/\t/,/g' > ${READ_SUFFIX}_heatmap_table_genus_no_plasmid.csv") | sed ':a;N;$!ba;s/\n/ /g' > command.txt
parallel -j1 < command.txt
cp ${READ_SUFFIX}_heatmap_table_genus_no_plasmid.csv ${HEATMAP_DIR}

#done!
```
The heatmap file was then edited to remove rows where no ARG had more than 2% of its links linking to that host (added to row called "Other").
```
for file in [SAMPLE]_heatmap_table_genus_no_plasmid.csv ; do remove2=$(tail -n +2 $file | while read line ; do all_less=1 && for n in $(echo "$line" | cut -d"," -f2- | sed 's/,/ /g') ; do if (( $n >= 0.02 )) ; then all_less=0 ; fi ; done && if (( all_less )) ; then echo -e all_less'\t'$(echo "$line" | cut -d"," -f2- | sed 's/,/\\t/g') ; else echo $line ; fi ; done | grep "all_less" | cut -f2- | awk '{for (i=1;i<=NF;i++) sum[i]+=$i;}; END{for (i in sum) print sum[i];}' | paste -sd, ) && head -1 $file > ${file/.csv/_removed2.csv} && new=$(tail -n +2 $file | while read line ; do all_less=1 && for n in $(echo "$line" | cut -d"," -f2- | sed 's/,/ /g') ; do if (( $n >= 0.02 )) ; then all_less=0 ; fi ; done && if (( all_less )) ; then echo -e all_less'\t'$(echo "$line" | cut -d"," -f2- | sed 's/,/\\t/g') ; else echo $line ; fi ; done) && echo "${new}" | grep -v "all_less" >> ${file/.csv/_removed2.csv} && echo "Other",$remove2 >> ${file/.csv/_removed2.csv} ; done
```
## Making the heatmaps
A heatmap of the ARG-host links for each sample was made using the [pheatmap](https://github.com/raivokolde/pheatmap) package in [R](https://www.r-project.org/).

Make the heatmap in R:
```
#install packages:

install.packages("pheatmap")
install.packages("RColorBrewer")

#_________________________________________________________________

#load packages:

library(pheatmap)
library(RColorBrewer)

#_________________________________________________________________

#select your colours (I recommend using 3, change colours by changing hex colour codes):

custom_colours <- c(
  "#0a2176", "#FFFFFF", "#FF7C12"
)

#get gradient for those colours:

breaksList = seq(0, 1, by = 0.01)
customcolour_gradient <- colorRampPalette(custom_colours)(length(breaksList))

#_________________________________________________________________

#load heatmap file:

[SAMPLE]_heatmap_table_genus_no_plasmid_removed2 <- read.csv(file="[SAMPLE]_heatmap_table_genus_no_plasmid_removed2.csv",header=TRUE,row.names = 1)

#OPTIONAL - if you want proper ARG name formats, create list with correct format for gene names e.g.:
#MUST BE IN SAME ORDER AS IN [SAMPLE]_heatmap_table_genus_no_plasmid_removed2

collist_args <- c(expression(paste(italic("tet"), "(M)_12")),
                       expression(paste(italic("aph(3')-III"), "_1")),
                       expression(paste(italic("tet"), "(32)_2")),
                       expression(paste(italic("bla"),""[OXA-500], "_1")),
                       expression(paste(italic("erm"), "(B)_12")),
                       expression(paste(italic("bla"),""[TEM-1],""[B], "_1")),
                       expression(paste(italic("msr"), "(D)_2"))
)

#same for host names, also MUST BE IN SAME ORDER AS IN [SAMPLE]_heatmap_table_genus_no_plasmid_removed2

rowlist_names <- c("Acutalibacteraceae",
                       bquote(italic("Agathobacter")),
                       bquote(italic("Eubacterium")),
                       bquote(italic("Faecousia")),
                       bquote(italic("Parabacteroides")),
                       "Other"
)

#now generate the heatmap (if you don't use the custom column/row names then add # before labels_row and labels_col):

pheatmap(
  method = c("pearson"),
  clustering_method = "complete",
  mat = [SAMPLE]_heatmap_table_genus_no_plasmid_removed2,
  breaks = breaksList,
# legend_breaks = c(0,1),
  treeheight_row = 50,
  treeheight_col = 50,
# fontsize_col = 10,
# fontsize_row = 10,
  fontsize = 8,
  cellwidth = 10,
  cellheight = 10,
  cluster_row = F,
  cluster_col = T,
  show_rownames = T,
  show_colnames = T,
  labels_row = as.expression(rowlist_names),
  labels_col = as.expression(collist_args),
# gaps_row = 1,
  legend = T,
  angle_col = 90,
  border_color = c("#0d0d0d"),
  color = customcolour_gradient
)
```

# Updated workflow
The workflow presented here was used in the meta3C/Hi-C manuscript. An updated workflow called [H-LARGe (Host-Linkage to antimicrobial resistance genes)](https://github.com/gregmcc97/H-LARGe) involving binning of the contigs for better host classification can be found [here](https://github.com/gregmcc97/H-LARGe).
