#!/bin/bash

#Trimming Single End Reads
SE_SAMPLES="H3K4me1_GM12878_input_rep1 H3K4me1_GM12878_input_rep2 H3K4me1_GM12878_rep1 H3K4me1_GM12878_rep2"

for SAMPLE in $SE_SAMPLES; do
    echo "${SAMPLE} trimming started"
    trim_galore -q 20 --stringency 2 --cores 4 -o ./qc/trim  --fastqc_args "--outdir ./qc/qc_trim --threads 4" ./data/reads/${SAMPLE}.fastq 
    echo "$SAMPLE trimming successful"

done


#Trimming Paired End Reads
PE_SAMPLES_1=(Nfxl1_GM12878_input_rep1_1 Nfxl1_GM12878_input_rep2_1 Nfxl1_GM12878_rep1_1 Nfxl1_GM12878_rep2_1 Nfxl1_K562_input_rep1_1 Nfxl1_K562_input_rep2_1 Nfxl1_K562_rep1_1 Nfxl1_K562_rep2_1)
PE_SAMPLES_2=(Nfxl1_GM12878_input_rep1_2 Nfxl1_GM12878_input_rep2_2 Nfxl1_GM12878_rep1_2 Nfxl1_GM12878_rep2_2 Nfxl1_K562_input_rep1_2 Nfxl1_K562_input_rep2_2 Nfxl1_K562_rep1_2 Nfxl1_K562_rep2_2)


for i in {0..7}; do
    echo "${PE_SAMPLES_1[i]} and ${PE_SAMPLES_2[i]} trimming started."
	trim_galore --paired -q 20 --stringency 2 --cores 4 -o ./qc/trim --fastqc_args "--outdir ./qc/qc_trim -threads 4" ./data/reads/${PE_SAMPLES_1[i]}.fastq ./data/reads/${PE_SAMPLES_2[i]}.fastq

done