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

