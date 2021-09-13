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

