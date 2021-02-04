#!/bin/bash

# Map RNAseq reads to the human transcriptome using Salmon

for i in `cat all.txt`;
do
	mkdir "./${i}" 

echo "Processing sample ${i}"
salmon quant -p 10 --gcBias --useVBOpt --seqBias --validateMappings -l A -i /media/NAS/nas-hug/ngs_data/cdna_homo_sapiens/release-101/Homo_sapiens.GRCh38.cdna.all.fa_index \
	 -1 /media/NAS/nas-hug/ongoing_projects/tnbc_unige_Intidhar_19-20/raw_data/${i}_R1.fastq.gz \
         -2 /media/NAS/nas-hug/ongoing_projects/tnbc_unige_Intidhar_19-20/raw_data/${i}_R2.fastq.gz \
         -o ./${i} \
	 -g /media/NAS/nas-hug/ngs_data/cdna_homo_sapiens/release-101/Homo_sapiens.GRCh38.101.gtf 
done 
