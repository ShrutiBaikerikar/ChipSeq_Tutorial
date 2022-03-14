#!/bin/bash

#Alignment of Single End Reads
SE_SAMPLES="H3K4me1_GM12878_input_rep1 H3K4me1_GM12878_input_rep2 H3K4me1_GM12878_rep1 H3K4me1_GM12878_rep2"

for SAMPLE in $SE_SAMPLES; do
    echo "${SAMPLE} alignment started"
    bowtie2 -p 4 -q -x ./align/index/GRCh38_noalt_as -U ./qc/trim/${SAMPLE}_trimmed.fq -S ./alignment/output/${SAMPLE}_unsorted.sam
    echo "$SAMPLE alignment successful"

done


#Alignment of Paired End Reads
PE_SAMPLES_1=(Nfxl1_GM12878_input_rep1_1_val_1 Nfxl1_GM12878_input_rep2_1_val_1 Nfxl1_GM12878_rep1_1_val_1 Nfxl1_GM12878_rep2_1_val_1 Nfxl1_K562_input_rep1_1_val_1 Nfxl1_K562_input_rep2_1_val_1 Nfxl1_K562_rep1_1_val_1 Nfxl1_K562_rep2_1_val_1)
PE_SAMPLES_2=(Nfxl1_GM12878_input_rep1_2_val_2 Nfxl1_GM12878_input_rep2_2_val_2 Nfxl1_GM12878_rep1_2_val_2 Nfxl1_GM12878_rep2_2_val_2 Nfxl1_K562_input_rep1_2_val_2 Nfxl1_K562_input_rep2_2_val_2 Nfxl1_K562_rep1_2_val_2 Nfxl1_K562_rep2_2_val_2)
PE_SAMPLES_3=(Nfxl1_GM12878_input_rep1 Nfxl1_GM12878_input_rep2 Nfxl1_GM12878_rep1 Nfxl1_GM12878_rep2 Nfxl1_K562_input_rep1 Nfxl1_K562_input_rep2 Nfxl1_K562_rep1 Nfxl1_K562_rep2)

for i in {0..3}; do
    echo "${PE_SAMPLES_1[i]} and ${PE_SAMPLES_2[i]} alignment started."
	bowtie2 -p 4 -q -x ./align/index/GRCh38_noalt_as -1 ./qc/trim/${PE_SAMPLES_1[i]}.fq -2 ./qc/trim/${PE_SAMPLES_2[i]}.fq -S ./alignment/output/${PE_SAMPLES_3[i]}_unsorted.sam
    echo "${PE_SAMPLES_1[i]} and ${PE_SAMPLES_2[i]} alignment successful."

done