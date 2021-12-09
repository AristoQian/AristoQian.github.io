#!/bin/bash
#SBATCH -A e31265
#SBATCH -p normal
#SBATCH -N 1
#SBATCH -n 12
#SBATCH -t 06:00:00
#SBATCH --mem=80gb
#SBATCH --job -name="opt5_chenxing"

cd $SLURM_SUBMIT_DIR
module purge all
module load STAR/2.7.9a

STAR --runThreadN 20 --soloType CB_UMI_Simple --soloCBwhitelist ./737K-august-2016.txt --genomeDir /projects/e31265/humanref --outFilterType BySJout --alignIntronMax 100000 --quantMode GeneCounts --outSAMtype None  --soloCellFilter EmptyDrops_CR --soloFeatures Gene --readFilesPrefix /projects/e31265/marrow/ --readFilesCommand "gzip - cd" --readFilesIn marrow_S1_L001_R2_001.fastq.gz,marrow_S1_L002_R2_001.fastq.gz,marrow_S1_L003_R2_001.fastq.gz,marrow_S1_L004_R2_001.fastq.gz,marrow_S1_L005_R2_001.fastq.gz,marrow_S1_L006_R2_001.fastq.gz,marrow_S1_L007_R2_001.fastq.gz,marrow_S1_L008_R2_001.fastq.gz marrow_S1_L001_R1_001.fastq.gz,marrow_S1_L002_R1_001.fastq.gz,
marrow_S1_L003_R1_001.fastq.gz,marrow_S1_L004_R1_001.fastq.gz,marrow_S1_L005_R1_001.fastq.gz,marrow_S1_L006_R1_001.fastq.gz,marrow_S1_L007_R1_001.fastq.gz,marrow_S1_L008_R1_001.fastq.gz