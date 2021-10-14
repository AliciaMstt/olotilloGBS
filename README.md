# Evolutionary plant breeding of Olotillo: Single Nuclueotide Variants (SNVs) calling

This repository contains the scripts to perform the SNV calling for the first 95 samples of Olotillo. This pipeline was standardized by Idalia Rojas based on the pipeline used for  the paper: [Contemporary evolution of maize landraces and their wild relatives influenced by gene flow with modern maize varieties](https://doi.org/10.1073/pnas.1817664116)

# Olotillo SNV Calling
##  Workflow and software requirements

1. Demultiplex: [GBSx](https://github.com/GenomicsCoreLeuven/GBSX). Local installation.

2. Alignment: Reads mapping against reference genome: [Picard tools](https://broadinstitute.github.io/picard/) run as a local java -jar installation. [Samtools](http://www.htslib.org/download/) and [NextGenMap](https://cibiv.github.io/NextGenMap/) are installed on Conabio's server.

3. Variant discovery: [GATK 4.2.1](https://github.com/broadinstitute/gatk/releases). Local installation.
 
4. Genotypes merging:  [GATK 4.2.1.0](https://github.com/broadinstitute/gatk/releases). Local installation.

## 01 Demultiplex

This step sort and separate the illumina reads  for multiple genotypes that were sequenced in the  same sequencing lane.

This step is executed with GBS and uses as input:
a) Multiples fastq file
b) barcode_file.txt: A tab separated text file whithout header,first column is the sample name, secon column is the barcode sequence and the third colum with the restriction enzyme i.e.:

```
GRO_N14_1	CTCC	ApeKI
CHIS_L36_1	TTCTC	ApeKI
CHIS_L43_4	GCTTA	ApeKI
CHIS_L42_1	AACGCCT	ApeKI
GRO_N13_2	AGGC	ApeKI
NAY_N27_3	TCGTT	ApeKI
ROO_E12_1	TGGCTA	ApeKI
GRO_N23_1	TGCTGGA	ApeKI
CHIS_L35_2	TGCA	ApeKI
GRO_N26_4	AGCCC	ApeKI

```

script: `01_demultiplex.sh`
For further  options check the GBSx github [repository](https://github.com/GenomicsCoreLeuven/GBSX)

```{bash}
#!/bin/bash
#SBATCH -p cluster
#SBATCH -n 1
#SBATCH --mem=20000
#SBATCH --job-name=stp_01_demult #Give your job a name.
#SBATCH --cpus-per-task=3 #Multithreading.
#SBATCH --error=job.%J.err #Std Error write standard error to this file
#SBATCH --output=job.%J.out #Std Out write standard output to this file
#SBATCH --mail-type=FAIL #Notify user by email when certain event types occur (BEGIN, END, FAIL, REQUEUE)
#SBATCH --mail-user=icrojasb@iecologia.unam.mx #Email for notifications from previous line

#This command was run to demultiplex samples from fastq files with GBSX_v1.3 
#12.09.21

#paths
GBSX="/LUSTRE/Genetica/common/olotillo/bin/GBSX"
bin="/LUSTRE/Genetica/common/olotillo/olotilloGBS/bin"
meta="/LUSTRE/Genetica/common/olotillo/olotilloGBS/meta"
data="/LUSTRE/Genetica/common/olotillo/olotilloGBS/genetic/"
out="/LUSTRE/Genetica/common/olotillo/olotilloGBS/out" #Just for the VCF file

#Note: Verify the barcode file does not contain and extra tab or enter

java -jar $GBSX/releases/latest/GBSX_v1.3.jar --Demultiplexer -f1 $data/plate1/Plate-1-20210721_S80_L003_R1_001.fastq.gz -f2 $data/plate1/Plate-1-20210721_S80_L003_R2_001.fastq.gz -i $meta/plate1_meta_demultiplex.txt -gzip true -o $data/fastq_demultiplexed/

```
## 02. Alignment and sorting

The demultiplexed reads were aligned against the B73 assembly version:  [Zm-B73-REFERENCE-NAM-5.0](https://download.maizegdb.org/Zm-B73-REFERENCE-NAM-5.0), released in January 2020. Check the MaizeGBD website for further details. 

script: `02_alignment.sh`
```{}
#!/bin/bash
#SBATCH -p cluster
#SBATCH -w nodo7
#SBATCH --mem=20000
#SBATCH --error=job.%J.err #Std Error write standard error to this file
#SBATCH --output=job.%J.out #Std Out write standard output to this file
#SBATCH --mail-type=FAIL #Notify user by email when certain event types occur (BEGIN, END, FAIL, REQUEUE)
#SBATCH --mail-user=icrojasb@iecologia.unam.mx #Email for notifications from previous line

#paths
home="/LUSTRE/Genetica/common/olotillo/olotilloGBS/"
bin="/LUSTRE/Genetica/common/olotillo/olotilloGBS/bin"
meta="/LUSTRE/Genetica/common/olotillo/olotilloGBS/meta"
data="/LUSTRE/Genetica/common/olotillo/olotilloGBS/genetic/"
out="/LUSTRE/Genetica/common/olotillo/olotilloGBS/out" #Just for the VCF file
fastq="/LUSTRE/Genetica/common/olotillo/olotilloGBS/genetic/fastq_demultiplexed"
ref="/LUSTRE/Genetica/common/olotillo/olotilloGBS/genetic/ref/Zm-B73-REFERENCE-NAM-5.0.fa"
bam="/LUSTRE/Genetica/common/olotillo/olotilloGBS/genetic/bam"
log="/LUSTRE/Genetica/common/olotillo/olotilloGBS/genetic/log"
project="olotillo-EvolBreed"
#ngm="/LUSTRE/Genetica/common/olotillo/bin/NextGenMap-0.5.0/build" NextGen Map is installed on CONABIO'S cluster
picardtools="/LUSTRE/Genetica/common/olotillo/bin"

#Create reference index
samtools faidx $ref 

#Create a list of files
ls $fastq |  grep fastq | egrep 'R1|R2'|sed s/.R1.fastq.gz// |sed s/.R2.fastq.gz// | sort -r |uniq > $meta/samplelist_$project.txt

while read name
do

    echo "Aligning paired $name"
    ngm -r $ref -1 $fastq/${name}.R1.fastq.gz -2 $fastq/${name}.R2.fastq.gz -o $bam/$name.paired.bam -t $ncores -b --rg-id $name --rg-sm $name --rg-pl illumina --rg-pu $project 
    
    echo "Aligning unpaired $name"
    ngm -r $ref -q $fastq/${name}.R1.fastq.gz -o $bam/$name.R1.unpaired.bam -t $ncores -b
    ngm -r $ref -q $fastq/${name}.R2.fastq.gz -o $bam/$name.R2.unpaired.bam -t $ncores -b
#
    echo "Processing for $name"

java -jar $picardtools/picard.jar MergeSamFiles I=$bam/$name.paired.bam I=$bam/$name.R1.unpaired.bam I=$bam/$name.R2.unpaired.bam O=$bam/$name.merged.bam VALIDATION_STRINGENCY=LENIENT COMPRESSION_LEVEL=0 
	
samtools view -h $bam/$name.merged.bam > $bam/$name.merged.sam
samtools view -hSb $bam/$name.merged.sam > $bam/$name.merged.v0.bam
        
java -jar $picardtools/picard.jar SortSam I=$bam/$name.merged.v0.bam O=$bam/$name.merged.v1.bam SORT_ORDER=coordinate VALIDATION_STRINGENCY=LENIENT COMPRESSION_LEVEL=0
java -jar $picardtools/picard.jar AddOrReplaceReadGroups I=$bam/$name.merged.v1.bam O=$bam/$name.merged.v2.bam \
RGID=$name RGLB=$project RGPL=ILLUMINA RGPU=$project RGSM=$name SORT_ORDER=coordinate VALIDATION_STRINGENCY=LENIENT COMPRESSION_LEVEL=0
java -jar $picardtools/picard.jar CleanSam I=$bam/$name.merged.v2.bam O=$bam/$name.merged.v3.bam VALIDATION_STRINGENCY=LENIENT
java -jar $picardtools/picard.jar MarkDuplicates I=$bam/$name.merged.v3.bam O=$bam/$name.bam  M=$bam/$name.marked_dup_metrics.txt
	rm $bam/$name.paired.bam
	rm $bam/$name.R1.unpaired.bam
	rm $bam/$name.R2.unpaired.bam
	rm $bam/$name.merged.sam
	rm $bam/$name.merged.bam
	rm $bam/$name.merged.v0.bam
	rm $bam/$name.merged.v1.bam
	rm $bam/$name.merged.v2.bam
	rm $bam/$name.merged.v3.bam

done < $meta/samplelist_$project.txt

```
## 03. Variant calling
We used GATK v.4.2.1.0 to perform the genetic varian calling and mark duplicates. This step  was run in parallel for each genotype.

The script `make_03_variantDiscovery.sh` receives as input a text file with a fastqc list, and create a bash script for each of them.
```{bash variant-calling}
while read name
do
echo -e '#!/bin/bash
 
#SBATCH -p cluster
#SBATCH -w nodo7
#SBATCH --mem=20000
#SBATCH --error=job.%J.err #Std Error write standard error to this file
#SBATCH --output=job.%J.out #Std Out write standard output to this file
#SBATCH --mail-type=FAIL #Notify user by email when certain event types occur (BEGIN, END, FAIL, REQUEUE)
#SBATCH --mail-user=icrojasb@iecologia.unam.mx #Email for notifications from previous line

#The Genome Analysis Toolkit (GATK) v4.2.1.0


home="/LUSTRE/Genetica/common/olotillo/olotilloGBS/"
bin="/LUSTRE/Genetica/common/olotillo/olotilloGBS/bin"
meta="/LUSTRE/Genetica/common/olotillo/olotilloGBS/meta"
data="/LUSTRE/Genetica/common/olotillo/olotilloGBS/genetic/"
out="/LUSTRE/Genetica/common/olotillo/olotilloGBS/out" #Just for the VCF file
fastq="/LUSTRE/Genetica/common/olotillo/olotilloGBS/genetic/fastq_demultiplexed"
bam="/LUSTRE/Genetica/common/olotillo/olotilloGBS/genetic/bam"
gvcf="/LUSTRE/Genetica/common/olotillo/olotilloGBS/genetic/gvcf"
ref="/LUSTRE/Genetica/common/olotillo/olotilloGBS/genetic/ref/Zm-B73-REFERENCE-NAM-5.0.fa"
log="/LUSTRE/Genetica/common/olotillo/olotilloGBS/genetic/log"
project="olotillo-EvolBreed"
gatk="/LUSTRE/Genetica/common/olotillo/bin/gatk-4.2.1.0/gatk"



#Create a dictionary for the reference genome
samtools dict Zm-B73-REFERENCE-NAM-5.0.fa > Zm-B73-REFERENCE-NAM-5.0.dict

#Create an index for each bam file
samtools index -b $bam/'$name'.bam \

#Call GATK HaplotypeCaller
$gatk --java-options "-Xmx20g" HaplotypeCaller  \
	-R $ref \
	-I $bam/'$name'.bam \
	-O $gvcf/'$name'.Olotillo.gvcf.vcf \
	-ploidy 2 \
	--do-not-run-physical-phasing \
	-ERC GVCF  #Emmitting reference confidence scores, gvcf format
	' > /LUSTRE/Genetica/common/olotillo/olotilloGBS/bin/03_variantDiscovery_$name.parallel.sh

done < /LUSTRE/Genetica/common/olotillo/olotilloGBS/meta/samplelist_olotillo-EvolBreed.txt

```

## 04.  CombineGVCFs and  GenotypeGVCFs

This step Combine per-sample gVCF files produced by HaplotypeCaller into a multi-sample gVCF file, check [GATK](https://gatk.broadinstitute.org/hc/en-us/articles/360037053272-CombineGVCFs) website for further details.

This is the most time-consuming step, for big genomes as maize run this per chromosome.

> Written with [StackEdit](https://stackedit.io/).
