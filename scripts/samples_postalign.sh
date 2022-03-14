#!/bin/bash

#SAM/BAM Processing of Single End Reads
SE_SAMPLES="H3K4me1_GM12878_input_rep1 H3K4me1_GM12878_input_rep2 H3K4me1_GM12878_rep1 H3K4me1_GM12878_rep2"

for SAMPLE in $SE_SAMPLES; do
    echo "${SAMPLE} post-alignment processing started"
    samtools view -h -S -b ./alignment/output/${SAMPLE}_unsorted.sam > ./alignment/unsorted_bam/${SAMPLE}_unsorted.bam
	sambamba-0.8.0 sort -t 4 -o ./alignment/sorted_bam/${SAMPLE}_sorted.bam ./alignment/unsorted_bam/${SAMPLE}_unsorted.bam
	sambamba-0.8.0 view -h -t 4 -f bam -F "[XS] == null and not unmapped and not duplicate" ./alignment/sorted_bam/${SAMPLE}_sorted.bam > ./alignment/filtered_bam/${SAMPLE}_aln.bam
    sambamba-0.8.0 sort -t 4 -o ./alignment/sorted_filtered_bam/${SAMPLE}_aln_sorted.bam ./alignment/filtered_bam/${SAMPLE}_aln.bam
	echo "$SAMPLE post-alignment processing successful"

done


#SAM/BAM Processing of Paired End Reads
PE_SAMPLES=(Nfxl1_GM12878_input_rep1 Nfxl1_GM12878_input_rep2 Nfxl1_GM12878_rep1 Nfxl1_GM12878_rep2 Nfxl1_K562_input_rep1 Nfxl1_K562_input_rep2 Nfxl1_K562_rep1 Nfxl1_K562_rep2)

for SAMPLE in $PE_SAMPLES; do
    echo "${SAMPLE} post-alignment processing started"
    samtools view -h -S -b ./alignment/output/${SAMPLE}_unsorted.sam > ./alignment/unsorted_bam/${SAMPLE}_unsorted.bam
	sambamba-0.8.0 sort -t 4 -o ./alignment/sorted_bam/${SAMPLE}_sorted.bam ./alignment/unsorted_bam/${SAMPLE}_unsorted.bam
	sambamba-0.8.0 view -h -t 4 -f bam -F "[XS] == null and proper_pair and not (unmapped or mate_is_unmapped) and not duplicate" ./alignment/sorted_bam/${SAMPLE}_sorted.bam > ./alignment/filtered_bam/${SAMPLE}_aln.bam
    sambamba-0.8.0 sort -t 4 -o ./alignment/sorted_filtered_bam/${SAMPLE}_aln_sorted.bam ./alignment/filtered_bam/${SAMPLE}_aln.bam
	echo "$SAMPLE post-alignment processing successful"

done