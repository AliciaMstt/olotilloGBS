#!/bin/bash
#SBATCH -p cluster
#SBATCH -w nodo7
#SBATCH --mem=60000
#SBATCH --error=job.%J.err #Std Error write standard error to this file
#SBATCH --output=job.%J.out #Std Out write standard output to this file
#SBATCH --mail-type=FAIL #Notify user by email when certain event types occur (BEGIN, END, FAIL, REQUEUE)
#SBATCH --mail-user=icrojasb@iecologia.unam.mx #Email for notifications from previous line

home="/LUSTRE/Genetica/common/olotillo/olotilloGBS/"
bin="/LUSTRE/Genetica/common/olotillo/olotilloGBS/bin"
meta="/LUSTRE/Genetica/common/olotillo/olotilloGBS/meta"
data="/LUSTRE/Genetica/common/olotillo/olotilloGBS/genetic/"
out="/LUSTRE/Genetica/common/olotillo/olotilloGBS/out60g" #Just for the VCF file
combinedGVCFs="/LUSTRE/Genetica/common/olotillo/olotilloGBS/genetic/gvcf"
ref="/LUSTRE/Genetica/common/olotillo/olotilloGBS/genetic/ref/Zm-B73-REFERENCE-NAM-5.0.fa"
log="/LUSTRE/Genetica/common/olotillo/olotilloGBS/genetic/log"
gatk="/LUSTRE/Genetica/common/olotillo/bin/gatk-4.2.1.0/gatk"
project="olotillo-EvolBreed"
tmp="/LUSTRE/Genetica/common/olotillo/olotilloGBS/tmp"


ls $combinedGVCFs | grep "vcf" | grep -v ".idx"   > $meta/mergeGVCFs.samplelist_$project.txt
tmp=""

while read prefix 
do
        tmp="$tmp --variant $combinedGVCFs/$prefix"
done < $meta/mergeGVCFs.samplelist_$project.txt

#For GATK 4
#CombineGVCFs is meant to be used for merging of GVCFs that will eventually be input into GenotypeGVCFs.
$gatk --java-options "-Xmx60g" CombineGVCFs \
   -R $ref \
   $tmp \
   -O $out/${project}_gatk4.1.9.comGty.vcf.gz \
   --read-filter MappingQualityReadFilter \
   --read-filter OverclippedReadFilter \
   
#GenotypeGVCF: perform joint genotyping on a single input

$gatk --java-options "-Xmx60g" GenotypeGVCFs \
 -R $ref \
 -V $out/${project}_gatk4.2.1.comGty.vcf.gz  \
 -ploidy 2 \
 -O $out/${project}_gatk4.2.1.vcf.gz
 
# rm $out/${project}_gatk4.2.1.comGty.vcf.gz 