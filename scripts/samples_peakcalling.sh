#!/bin/bash

#Peak Calling For Single End Reads
SE_SAMPLES_1="H3K4me1_GM12878_rep1 H3K4me1_GM12878_rep2"
SE_SAMPLES_2="H3K4me1_GM12878_input_rep1 H3K4me1_GM12878_input_rep2"

for j in {0..1}; do
    echo "${SE_SAMPLES_1[j]} and ${SE_SAMPLES_2[j]} peak calling started."
	macs2 callpeak -t ./alignment/filtered_bam/${SE_SAMPLES_1[j]}_aln.bam -c ./alignment/filtered_bam/${SE_SAMPLES_2[j]}_aln.bam -f BAM --broad -g hs --broad-cutoff 0.1 --keep-dup auto --outdir ./macs2/output -n ${SE_SAMPLES_1[j]}
    echo "${SE_SAMPLES_1[j]} and ${SE_SAMPLES_2[j]} peak calling successful."

done


#Peak Calling For Paired End Reads
PE_SAMPLES_1=(Nfxl1_GM12878_rep2 Nfxl1_K562_rep1 Nfxl1_K562_rep2)
PE_SAMPLES_2=(Nfxl1_GM12878_input_rep2 Nfxl1_K562_input_rep1 Nfxl1_K562_input_rep2)

for i in {0..2}; do
    echo "${PE_SAMPLES_1[i]} and ${PE_SAMPLES_2[i]} peak calling started."
	macs2 callpeak -t ./alignment/filtered_bam/${PE_SAMPLES_1[i]}_aln.bam -c ./alignment/filtered_bam/${PE_SAMPLES_2[i]}_aln.bam -f BAMPE -g hs --keep-dup auto --outdir ./macs2/output -n ${PE_SAMPLES_1[i]} -B
	echo "${PE_SAMPLES_1[i]} and ${PE_SAMPLES_2[i]} peak calling successful."

done