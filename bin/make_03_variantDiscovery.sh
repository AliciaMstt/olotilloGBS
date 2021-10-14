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

