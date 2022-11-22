# ChIP-Sequencing Analysis of H3K4me1 and Nfxl1: Motif Analysis Tutorial

Next Generation Sequencing comprise several technologies that help us understand genome wide expression data and provide high resolution insights even up to the transcript activity. ChIP sequencing or Chromatin Immunoprecipitation followed by high-thoroughput sequencing is one such technique that focuses on protein-DNA interactions.

It aids in identification of genomic sites to which the proteins of interest bind. These proteins generally include histone marks or transcription factors which regulate transcription and hence decide which genes will be expressed and which won't.

Differential Binding and Motif Analysis are the common techniques used in downstream analysis of ChIP-Seq data. It helps us identify which are the different binding sites to which the transcription facor or histone mark binds under different conditions (i.e. disease/treatment vs reference/normal conditions) and in turn the different genes whose expression it may influence. 

Motif analysis helps identify the underlying sequence of the binding site.

In this tutorial, we will be going over a general workflow for Motif Analysis using ChIP Sequencing data. While the procedures demonstrated in this tutorial are the most common steps involved in ChIP Sequencing data analysis, there are a few things you must keep in mind:
* There are many tools available for each step of ChIP Sequencing analysis. You can choose any of them and this choice is heavily influenced by type of data, organism, sample number, statiscal requirements as well as goals of research/analysis.
* ChIP Sequencing data analysis is not limited to identifying differential binding sites alone; it has many applications and these require separate workflows.
* The tools used in this tutorial help analyse the chosen dataset but this may not be the same for you. You may have to choose separate tools for your ChIP-Seq data.
* When working with large number of samples, it is advisable to work on a high performance cluster or utilise cloud services.

Let's begin the tutorial with a basic introduction to ChIP Sequencing.

## Table of Contents

- [Introduction](#intro)
     - [Introduction to ChIP sequencing](#general_intro)
     - [Background of the dataset](#dataset_intro)
     - [Installation](#installation)
     - [Setting up Working Directory](#work_directory)
- [Motif Analysis from ChIP-Seq Data](#workflow)
     - [Quality Control of ChIP-Seq Reads](#quality_control)
     - [Preprocessing of ChIP-Seq Reads](#preprocessing)
     - [Alignment of ChIP-Seq Reads to Reference Genome](#alignment)
     - [SAM/BAM Processing of Aligned Reads](#sam_bam)
     - [Peak Calling of Aligned Reads](#peak_calling)
     - [ChIP specific quality control](#chip_qc)
     - [Handling replicates with IDR and bedTools](#replicates)
     - [Visualization](#visualization)
     - [Functional Analysis](#functional)
     - [Motif Analysis](#motif)
- [Conclusion](#conclusion)
- [Citations](#citations_list)
- [License](#license_name)

## Introduction <a name="intro"></a>

### What is ChIP Sequencing? <a name="general_intro"></a>

ChIP sequencing or Chromatin Immunoprecipitation followed by high-throughput sequencing is a technique that identifies regions in the genome to which certain proteins bind; these proteins are Transcription Factors and post-translational histone modifications.

Transcription is a process that converts DNA into messenger RNA that can be further translated into protein. Alternative splicing during transcription determines which part of the gene will contribute to the final mRNA that will code for the protein variant. 

Protein-DNA interactions help regulate gene expression. Transcription can be regulated before initiation or even during the process to determine how many transcripts (mRNA) are produced or which variant of the transcript is produced.

The rate of transcription initiation is one of the primary ways by which gene expression is regulated. Transcription factors are proteins that bind to specific areas (motifs) in the DNA and regulate the transcription rate of that gene. They may enhance or lower transcription rate.

Regulation of DNA packing and structure can also affect gene expression. Changes in DNA packing in the nucleus can determine which regions of the DNA would be accessible to transcription factors and thus influence the rate of transcription.

Histone is a protein that helps in packing DNA in the nucleus. They serve as spools for thread-like DNA to bind around, thus condensing the structure into chromatin that can be packed in the nucleus.

Modifications in histone tails (such as acetylation, methylation and phosphorylation) can affect the affinity of histones to bind to DNA and thus impact DNA packing. Histone modifications are associated with a number of different transcription-related conditions.
They may limit or expose transcription factor binding sites thus indirectly regulating transcription and gene expression.

ChIP sequencing thus helps us study such protein-DNA interactions that can affect gene expression. 

The wet-lab protocol of ChIP sequencing begins with formaldehyde crosslinking of chromatin-bound proteins to DNA. This gives an idea of the protein and histone modifications that are distributed along the genome at that given time.

The crosslinked chromatin is sheared such that average size of the DNA fragment is around 200 bp. The next step involves immunoprecipitation. An antibody specific to the protein or modification of interest is introduced.

It binds to the protein-DNA complexes and these complexes are further enriched using beads that bind to the antibody.

The complexes are washed, proteins are digested and the DNA is extracted and further used in library preparation and high throughput sequencing.

After sequencing, we get the raw reads or data about the DNA sequence that is further analysed with bioinformatic and computational techniques.

In the next section, we will go over the different steps involved in the computational aspect of ChIP-Seq analysis. 


### What are the general steps involved in ChIP Sequencing Data Analysis?

Millions of reads are obtained from the ChIP sequencing experiments which are analysed with the help of computational approaches and statistics. These techniques assist in motif discovery, differential binding analysis, gene ontology and pathway analysis; all of which contribute to better understanding of transcriptional machinery, chromatin structure and gene expression regulation.

The following image gives you a brief idea of the computational approaches involved in ChIP-Seq data analysis.


<p align="center">
<img src="https://github.com/ShrutiBaikerikar/ChipSeq_Tutorial/blob/main/images/chip_seq_workflow.png" width="400" height=600 alt="ChIP Sequencing Workflow image"/>
</p>

<p align="center">
     <b>ChIP Sequencing Data Analysis - Workflow </b>
</p>

On the left are the basic steps involved in ChIP Sequencing Data Analysis and on the right are the tools/software that are used in these procedures. 


### What tools are will be using in this RNA Sequencing Analysis - Tutorial?

For this tutorial, we will be using the following tools
-	Quality Control : FASTQC, ChIPQC
-	Pre-processing: Trim Galore
-	Read Alignment: Bowtie2
-	SAM files Processing: Samtools, Sambamba
-	Peak Calling: Macs2
-	Motif Analysis: Homer

Please download these tools and install them in your system.
Before proceeding with the tutorial, please install Ubuntu if you are using Windows OS. We will be using Bash/Shell and RStudio in this tutorial.

### Which dataset will we use for DE gene analysis from RNA Seq data? <a name="dataset_intro"></a>

We will be using the following ENCODE datasets for our analysis:
* [ENCSR000AKF](https://www.encodeproject.org/experiments/ENCSR000AKF/) : H3K4me1 ChIP-Seq on GM12878
* [ENCSR746XEG](https://www.encodeproject.org/experiments/ENCSR746XEG/): NFXL1 ChIP-Seq on GM12878 
* [ENCSR085DD1](https://www.encodeproject.org/experiments/ENCSR085DDI/): NFXL1 ChIP-Seq on K562


H3K4me1 is a histone modification in the histone protein H3. A methyl group is added to the fourth lysine residue in H3. H3K4me1 enriched regions correspond to enhancers and promoters (sites on the DNA to which transcription factors bind and increase rate of transcription).

NFXL1 stands for Nuclear Transcription Factor, X box binding like 1. It enables DNA-binding transcription activity. Diseases associated with NFXL1 include speech and communication disorders.

H3K4me1 dataset has single end reads while both NFXL1 datasets have paired end reads.

You can download the respective datasets from ENCODE. However, it is important to note that while we will be using these datasets in the tutorial, the files are large and difficult to process on an average system or PC.

It would be advisable to go through the initial steps of the tutorial rather than implementing them for each file. You can directly utilise the BED files generated after peak calling and IDR and run the downstream analysis.

---------------------------------------------------------------------------------
### Installation of Packages <a name="installation"></a>

This tutorial assumes that you have the basic knowledge of Linux/Shell scripting and R programming. To implement this tutorial, please ensure that you have the following installed in your system:

- If you are using a Windows system, install Ubuntu
- R 4.0.0
- Python 3
- RStudio
- FastQC
- Bowtie2
- MACS2
- SAMtools, Sambamba
- Trim Galore
- Homer

The R packages that will be required later on in the analysis can be installed using RStudio and the installation procedure has been described further on in the tutorial.

--------------------------------------------------------------------------------------

### Setup your working directory <a name="work_directory"></a>

In the image above, it is pretty obvious that ChIP Seq data analysis consists of multiple steps involving different input and outputs. Therefore it is very important to organize your data into separate folders, so that you don't get lost while working. :)

Here is how you can organize your working directory structure:

``` bash
── chip_seq_analysis/
  │   └── alignment/                     <- Data generated during alignment steps
  │       ├── 1_index/                   <- Folder to store the indexed genome files from Bowtie2
  │       ├── 2_output/                  <- Alignment files generated from Bowtie2 (.SAM)
  │       ├── 3_unsorted_bam/            <- Aligned SAM Files converted to BAM (.BAM)
  │       ├── 3_sorted_bam/              <- Sorting and Indexing BAM Files (.BAM and .bai)
  │       ├── 3_filtered_bam/            <- Filtering BAM Files (.BAM)
  │       ├── 3_sorted_filtered_bam/     <- Resorting and indexing filtered BAM files (.BAM and .BAI)
  │       ├── 3_pooled_bam/              <- Generating pooled BAM files for replicate analysis (.BAM)
  │       
  │   
  │    └── bedtools/                     <- Data generated during replicate handling for Broad Peaks
  │   
  │    └── chipseq_r/                    <- Data generated during analysis of ChIP-Seq data in R
  │       ├── 1_bam/                     <- Folder to store the sorted and filtered .BAM files and their indexes (.BAM and .BAI)
  │       ├── 2_peaks/                   <- Folder to store peak files generated by MACS2 (.BED)
  │       ├── 3_consensus_peaks/         <- Folder to store consensus peak files generated by IDR and Bedtools (.BED)
  │
  │    └── data/                         <- Location of input files
  │       ├── 1_blacklist/               <- Folder to store the bed file for GRCh38 blacklisted regions (.BED)
  │       ├── 2_chr_bed/                 <- Folder to store the bed file which contains coordinates for all genes in the genome GRCh38
  │       ├── 3_reads/                   <- Folder to store all input and control ChIP-Seq data (.FASTQ)
  │       
  │   
  │    └── idr/                          <- Data generated during replicate handling for Narrow Peaks
  │         
  │    └── macs2/                        <- Data generated during peak calling
  │       ├── 1_output/                  <- Folder to store the bed files generated during peak calling under stringent conditions (.BED)
  │       ├── 2_output2/                 <- Folder to store the bed file generated during peak calling under relaxed conditions (.BED)
  │       
  │   
  │    └── motif/                        <- Data generated during motif analysis
  │       ├── 1_h3k4me1_gm12878_output/  <- Folder to store motif analysis output files for H3K4me1_GM12878 samples
  │       ├── 2_nfxl1_gm12878_output/    <- Folder to store motif analysis output files for Nfxl1_GM12878 samples
  │       ├── 3_nfxl1_k562_output/       <- Folder to store motif analysis output files for Nfxl1_K562 samples
  │  
  │    └── qc/                           <- Data generated during quality control and pre-processing steps
  │       ├── 1_fqc_results/             <- Results of FASTQC for each sample
  │       ├── 2_qc_trim/                 <- Results of FASTQC for every trimmed sample
  │       ├── 3_trim/                    <- Output files of Trimmed reads for every sample
  │
  │    └── visualization/                <- Data generated during visualization of ChIP-seq data
  │       ├── 1_bigwig/                  <- Folder to store BigWig files for each sample (.BW)
  │       ├── 2_figures/                 <- Folder to store visualization images for each sample
  │       ├── 3_matrix/                  <- Folder to store intermediate count matrix file for each sample
  │      
  │       
   
``` 
--------------------------------------------------------------------------------------------------

## Motif Analysis from ChIP Sequencing Data - Workflow <a name="workflow"></a>

Now, we will begin our data analysis. One important thing to note is that these FASTQ files contain a lot of data that can be difficult to process on average PC. Therefore it's best to use an HPC or cloud computing environment to run this tutorial.
I will soon be adding a toy dataset that can be used to implement this tutorial. Till then this tutorial best serves as a guide or reference when working with real data.

Lets's begin our analysis with quality control of our raw ChIP-Seq reads.

---------------------------------------------------

### 1. Quality Control of ChIP-Seq Samples <a name="quality_control"></a>

General quality control of high throughput sequencing reads is the first step in ChIP-Seq data analysis. It can help avoid problems that would occur during genome alignment. 
Apart from general quality control, ChIP-specific quality measures are conducted in the later stages of the analysis

Quality-related issues generally arise during sequencing or library preparation. These include: low confidence bases, PCR artifacts, sequence-specific bias, sequence contamination, untrimmed adapters, 3’/5’ positional bias.

FASTQC is a Java program that performs multiple quality checks on tens of millions of reads in a few minutes. It reports and helps visualize information on base content and quality, k-mer content, presence of ambiguous bases, overrepresented sequences, duplicates etc.

The following command produces a quality report for a single sample:

```bash

fastqc ./data/reads/H3K4me1_GM12878_input_rep1.fastq -o ./qc/fqc_results

```

You can run the same command (but replacing the name of the FASTQ file) for each sample or you can run this command for all samples at once

```bash

fastqc ./data/reads/*.fastq -o ./qc/fqc_results

```


Let’s look at some quality check reports produced by FASTQC for sample H3K4me1_GM12878_input_rep1.fastq
<p> </br> </p>
<p align="center">
<img src="https://github.com/ShrutiBaikerikar/ChipSeq_Tutorial/blob/main/images/H3K4me1_GM12878_input_rep1_fastqc_image_1.png" width="800" height=400 alt="FastQC basic statistics image"/>
</p>

<p align="center">
     <b>FastQC Report: Basic Statistics for sample H3K4me1_GM12878_input_rep1.fastq  </b>
</p>

At the left of the image, FASTQC has given judgements (pass, warn, fail) on several quality metrics. These judgements are based on general thresholds and poor judgements may not always suggest that the sample has failed quality checks.

In the Basic Statistics section, we can see that the sample H3K4me1_GM12878_input_rep1.fastq has 10676160 sequences with each read length of 51 bases. Also the base quality is given in Sanger/Illumina 1.9 encoding or Phred 33 score.

<p> </br> </p>
<p align="center">
<img src="https://github.com/ShrutiBaikerikar/ChipSeq_Tutorial/blob/main/images/H3K4me1_GM12878_input_rep1_fastqc_image_2.png" width="800" height=400 alt="FastQC per base sequence quality image"/>
</p>

<p align="center">
     <b>FastQC Report: Per base sequence quality for sample H3K4me1_GM12878_input_rep1.fastq </b>
</p>

This plot shows the range of quality values across all bases at each position in the FASTQ file. Base quality indicates the confidence in the base call or how correctly the sequencer has identified the base at a given position in the given sequence.

These scores are expressed in Phred scale which is given as Q = -10 log10 P where P is the probability that the base is wrong. The scores range from 0 to 40. In the FASTQ files, they are encoded as ASCII characters instead of numbers to save space.
Phred 33 score refers to the encoding where the 33rd ASCII character is used as 0. Phred 64 score is used in old Illumina software the 64th ASCII character is used as 0.

In the Base quality report, at each position a BoxWhisker type plot is drawn. Here the upper and lower whiskers indicate 10% and 90% points while blue line indicates mean quality and the red line indicates median quality.
The background of the graph divides the y axis into very good quality calls (green), calls of reasonable quality (orange), and calls of poor quality (red). 

For the sample H3K4me1_GM12878_input_rep1.fastq, we can see that quality scores are declining towards the 3’ end. This is a common phenomenon for Illumina sequencing reads and is attributed to different errors occurring in the sequencing process. 
This can be rectified by filtering reads with low average base quality or trimming low quality reads.

Other important metrics include:
* **Per base sequence composition plot**: It reports the percent of bases called at each position across all reads in the file. You can expect read start sequence biases in this plot if fragmenting with transposases or due to random hexamer priming.

* **Per Sequence GC content plot**:  The ‘Per Sequence GC content’ averages GC content over all sequences (indicated as red line) and compares it modelled normal distribution of GC content. 
In this sample, we observe that both the distributions are almost similar. In case, the red distribution would be unusually different from the blue one, this could mean that the library is contaminated with a genome of another organism (could be observed by broad peaks) or there are other kinds of bias (in case of over-represented sequences, you would see sharp peaks). 

* **Sequence Duplication Plot**: The ‘Sequence Duplication Levels’ plot shows the relative number of sequences with different degrees of duplication.
This module analyses only first 100,000 sequences in each file. Each sequence is tracked to the end of the file to give a representative count of the overall duplication level.
In a diverse library most sequences will occur only once in the final set. A low level of duplication may indicate a very high level of coverage of the target sequence, but a high level of duplication could indicate a bias such as PCR over-amplification or low complexity library due to small amounts of starting materials.
In this sample, we do not observe high level of duplication. Duplicates are removed prior to peak calling.

* **Over-represented sequences plot**: The ‘Overrepresented Sequence Plot’ lists all of the sequence which make up more than 0.1% of the total. 
Theoretically, a normal-high throughput library would have a diverse set of sequences; no individual sequence would account for a high fraction of the whole. However, if the sample does contain overrepresented sequences, it could mean that the plot is highly biologically significant or the library is contaminated or it has a bias.
With ChIP-Seq, it is quite likely to see over-represented sequences in the immunoprecipitated sample as you are enriching for specific protein-bound DNA sequences. Lack of over-represented sequences, as in case of this sample, doesn’t necessarily mean that you have a bad sample/experiment.

----------------------------------------------------------------

### 2. Preprocessing of ChIP-Seq samples <a name="preprocessing"></a>

After conducting quality checks, multiple pre-processing steps can be conducted to mitigate some of the quality problems that arose during the experimental setup. This assists in better alignment of the reads to the genome.
These steps include:
* **Filtering**: You can filter reads based on their quality. Average read quality is calculated and if it falls below the user-defined threshold that read is dropped from further analysis. In case of Paired End reads, if the mean quality of either of the reads drops below the given threshold, the read pair is dropped from further analysis.

* **Trimming**: Instead of dropping entire reads or read pairs from further analysis, you can trim low quality bases from a given read. There are several ways to trim a read. Some common techniques include trimming reads from either 3’ or 5’ end or using a sliding window approach where you can define a length of a search window and examine the mean quality of bases in that window. 

If the mean quality is above the given threshold, the window slides further to examine the next set of bases.

Sliding the window from the 5′ end keeps the beginning of the read until the quality falls below the defined threshold, while sliding from the 3′ end cuts until it reaches a window with good enough quality. The window size is an essential parameter that needs to be tuned; very small window size may lead to a stringent check and lead to loss of reads.

* **Removal of Adapters**: Adapters are short, known sequences of oligonucleotides that are used to extract or fish out DNA sequences of interest. Other tags such as primers and multiplexing identifiers may also be attached to the reads. These need to be removed prior to further analysis. However, removal of adapter content can present several hurdles.

Like any other part of the read, these tags can undergo sequencing errors like mismatches, indels and ambiguous bases. When sequencing small-RNAs, reads can run into a ‘read-through’ situation where the reads extend into the adapter and 3’ end adapter can be partial.

Trimming tools overcome these challenges by examining the reads for known adapter content, aligning the portions of the reads that partially match these adapter sequences and then trimming the matched portions.

* **Examining Read Length**: Read length distribution gives us an idea of how useful the reads are for further steps such as genome alignment, transcriptome assembly and detection of splice isoforms. Very short reads, resulting from trimming and adapter removal, map unambiguously to the genome and hence can be removed from further analysis.

The above-mentioned steps are the most common procedures utilised for pre-processing data prior to alignment. Depending on your data quality and the QC plots, you may need to take additional steps such as cautious removal of duplicates, examining reads for possible contaminants and removing related sequences etc.

Several tools are available for pre-processing the data: Trimmomatic, PRINSEQ, Cutadapt, TrimGalore, FastX, TagCleaner 

In this tutorial, we will be using Trim Galore. We will remove adapters and perform basic trimming and then examine the QC reports once again.

Assuming you have downloaded and installed Trim Galore. (https://www.bioinformatics.babraham.ac.uk/projects/trim_galore/) Add Trim Galore to your PATH or configuration file.

**For Single End Reads**

```bash

trim_galore -q 20 --stringency 2 --cores 4 -o ./qc/trim  --fastqc_args "--outdir ./qc/qc_trim --threads 4" ./data/H3K4me1_GM12878_input_rep1.fastq

```

**For Paired End Reads**

```bash

trim_galore --paired -q 20 --stringency 2 --cores 4 -o ./qc/trim --fastqc_args "--outdir ./qc/qc_trim -threads 4" ./data/reads/Nfxl1_GM12878_input_rep1_1.fastq ./data/reads/Nfxl1_GM12878_input_rep1_2.fastq

```

In this command, the options are:
- paired : It performs trimming for paired-end reads. After adapter removal and trimming, the read length is examined. If both reads from the paired group have length above the given threshold, they are retained.
If one only one read has length more than the given threshold then the read pair is discarded by default (this can be changed by specifying --retain_unpaired). This helps discard short reads.
The minimum length can be specified by --length parameter; by default the threshold is 20.
- q : It is the Phred quality score. Reads ends with quality score less than that specified with q are trimmed in addition to adapter removal. The default value is 20.
- stringency: It specifies the number of base pairs that need to overlap with the adapter content before it is trimmed from the sequence. When not specified, Trim Galore uses the first 13 bp of Illumina adapter 'AGATCGGAAGAGC'.
- cores: The number of threads available on your system for computing.
- o: It specifies the path to the output directory where the trimmed reads will be stored.
- fastqc_args: It specifies the extra arguments that re to be passed to FASTQC such as output directory and threads etc. Specifying --fastqc_args helps perform quality control check on trimmed reads.


You have to run the above mentioned commmand for each of your samples; just replace the name of the sample. Or you could run the bash script to process all samples together.
This script is available as samples_trim.sh in <u> <b> the scripts folder in the repository. </b> </u>

After trimming you will get files such as H3K4me1_input_rep1_trimmed.fq for single end reads and Nfxl1_GM12878_input_rep1_1_val_1.fq  and Nfxl1_GM12878_input_rep1_2_val_2.fq for paired end reads.

Open the FASTQC html file for each trimmed sample. We will check the results for H3K4me1_GM12878_input_rep1_trimmed.fq.

<p> </br> </p>
<p align="center">
<img src="https://github.com/ShrutiBaikerikar/ChipSeq_Tutorial/blob/main/images/H3K4me1_GM12878_input_rep1_trimmed_fastqc_image_3.png" width="800" height=400 alt="FastQC per base sequence quality for trimmed reads image"/>
</p>

<p align="center">
     <b>FastQC Report: Per base sequence quality for sample H3K4me1_GM12878_input_rep1_trimmed.fq </b>
</p>

After running FASTQC again on the trimmed reads, we observe slight improvement in our Sequence Quality per position plot. Also, there are changes in the ‘Sequence Length distribution’ plot because after trimming, reads are likely to be of different lengths.

You could consider trimming from the 5’ ends to further improve overall quality. However, we need to keep in mind that subsequently alignment of the reads to reference genome is to be performed using Bowtie2. And large proportion of reads shorter than 51 bp are unlikely to help in alignment.

You could either stick to this minimal trimming and adapter removal procedure performed above, or use soft clipping option during alignment with Bowtie2 or you could consider aligning the reads with BWA or Bowtie1.
For now, we will stick to minimal trimming and adapter removal.

-----------------------------------------------------------------

### 3. Genome Alignment of Preprocessed ChIP-Seq reads <a name="alignment"></a>

Genome alignment involves mapping reads to a reference genome. This helps us identify where the reads originate from and how similar they are to the reference genome.

However, there are a couple of challenges when mapping reads to a genome:
• There are millions of reads and they are short in length. Genomes are large and do contain sequences like repeats or pseudogenes which make it difficult to map each read to a unique position.
• Additionally, there are mismatches and indels caused by genomic variation and sequencing errors that need to be dealt with.

ChIP-Seq reads are derived from specific, short, enriched DNA fragments and thus the data is usually aligned with ‘short contiguous read mappers’.

The basic approach of such aligners is ‘seed and extend’ which involves:
1. identifying segments of reads of defined lengths (seeds) that precisely map to a given location in the genome. Seed matches can be exact or tolerate mismatches. Shorter seeds increase sensitivity while longer seeds permit faster searches.
2. extend the reads in both directions to map the rest of the read or maximum mappable length

Most aligners assign a quality score that estimates the accuracy of the obtained alignment.

One such aligner is Bowtie2. Some of the advantages Bowtie2 offers are that it is fast, has high accuracy ,and has low memory requirements.  
It is based on Burrows Wheeler Transform (BWT) which is a data compression algorithm. It transforms data in a way that it is amenable to compression and is useful for data containing lots of repeats (sequence information, in our case).
It is useful for aligning reads of about 50 up to 100s base pairs in length.

Bowtie2 supports gapped alignment with affine gap penalties. It offers end-to-end alignment where it searches for alignments including all characters in the given read as well as local alignment where it trims some characters from either or both ends to improve the alignment score.

You can download and install Bowtie2 from https://github.com/BenLangmead/bowtie2. 

Additionally, download the reference genome index (H. sapiens, GRCh38 no-alt analysis set) given by Bowtie2 developers to save time. You can download it from http://bowtie-bio.sourceforge.net/bowtie2/manual.shtml

Unzip and save the downloaded ‘genome’ index in the index folder under the alignment directory.

Add Bowtie2 to your path, if not added to the configuration file. Now run Bowtie2 using the following command:

**For Single End Reads**

```bash

bowtie2 -p 4 -q -x ./align/index/GRCh38_noalt_as -U ./qc/trim/H3K4me1_GM12878_input_rep1_trimmed.fq -S ./align/output/H3K4me1_GM12878_input_rep1_unsorted.sam

```

**For Paired End Reads**

```bash

bowtie2 -p 4 -q -x ./align/index/GRCh38_noalt_as -1 ./qc/trim/Nfxl1_GM12878_input_rep1_1_val_1.fq -2 ./qc/trim/Nfxl1_GM12878_input_rep1_1_val_1.fq -S ./align/output/Nfxl1_GM12878_input_rep1_unsorted.sam

```

In these commands, the parameters used are:
-p : The number of threads available for computing on your system
-q : Reads that are FASTQ files.
-x : The basename of the index for the reference genome.
-U: Files containing unpaired reads to be aligned
-1 and 2: Files containing mate 1s and 2s respectively for paired end reads
-S: File to write SAM alignments to.

SAM files or Sequence Alignment Map are large files and need to be processed further before peak calling. Let’s go over that in the next section.

You have to run this command for each sample by changing the sample names. Alternatively, you can use the samples_align.sh script present in **the scripts folder in the repository.** 

---------------------------------------------------------------------

### 4. SAM/BAM Processing of Aligned ChIP-Seq Reads <a name="sam_bam"></a>

The output from the Bowtie2 aligner is in the SAM format. SAM or Sequence Alignment Map format is a tab-delimited text file that stores information regarding the alignment of an individual read to the reference genome.

The file begins with the header section and each heading begins with the ‘@’ symbol. The header contains information about the source of data, reference genome used, method of alignment ,and other additional information about the alignment.

Following the header is the alignment section. It contains information for each read and where and how it aligns to the reference genome. The alignment section has 11 mandatory fields:
- QNAME: Query name or read name
- FLAG: numerical value that tells us about the mapping and whether the read is part of a pair
- RNAME: Reference sequence name
- POS: 1- based leftmost mapping position
- MAPQ: numerical value indicating alignment or mapping quality
- CIGAR: a sequence of letters and numbers that indicate which bases align to the reference (i.e. match, mismatch, deletion, insertion), and the numbers indicate the associated base lengths for each match/mismatch/insertion/deletion.
- RNEXT: Reference name of the mate/next read
- PNEXT: Position of the mate/next read
- TLEN: observed template length
- SEQ: segment sequence
- QUAL: ASCII of Phred-scaled base quality score for each base of the sequence

Apart from these mandatory fields, there are multiple optional fields or TAGS available that store aligner-specific information. 

Some attributes of the SAM format that are relevant to this tutorial are FLAG, MAPQ ,and among the optional fields [XS] Tag that is specifically used by Bowtie2.

The Bitwise FLAG helps us filter out reads that are uniquely mapped and concordantly mapped in the case of paired-end reads. Similarly, you can consider filtering even based on the Mapping Alignment Score (MAPQ); a threshold of 10-30 is generally used.

First, we need to convert SAM files to BAM files. BAM or Binary Alignment Map is the binary equivalent of SAM files and stores the same data in a compressed manner.
Next, we sort the BAM file by genomic coordinate before using them further. 
Finally, we filter out uniquely mapping reads. We discard multi-mapping reads to reduce false positives and increase confidence in site discovery.

We will use SAMtools and Sambamba for this step. You can download it from http://www.htslib.org/ and https://lomereiter.github.io/sambamba/ respectively.

Add these to your path before executing the following commands.

**For Single End reads**

```bash

samtools view -h -S -b ./alignment/output/H3K4me1_GM12878_input_rep1_unsorted.sam > ./alignment/unsorted_bam/H3K4me1_GM12878_input_rep1_unsorted.bam
sambamba-0.8.0 sort -t 4 -o ./alignment/sorted_bam/H3K4me1_GM12878_input_rep1_sorted.bam ./alignment/unsorted_bam/H3K4me1_GM12878_input_rep1_unsorted.bam
sambamba-0.8.0 view -h -t 4 -f bam -F "[XS] == null and not unmapped and not duplicate" ./alignment/sorted_bam/H3K4me1_GM12878_input_rep1_sorted.bam > ./alignment/filtered_bam/H3K4me1_GM12878_input_rep1_aln.bam
sambamba-0.8.0 sort -t 4 -o ./alignment/sorted_filtered_bam/H3K4me1_GM12878_input_rep1_aln_sorted.bam ./alignment/filtered_bam/H3K4me1_GM12878_input_rep1_aln.bam

```

**For Paired End reads**

```bash

samtools view -h -S -b ./alignment/output/Nfxl1_GM12878_input_rep1_unsorted.sam > ./alignment/unsorted_bam/Nfxl1_GM12878_input_rep1_unsorted.bam
sambamba-0.8.0 sort -t 4 -o ./alignment/sorted_bam/Nfxl1_GM12878_input_rep1_sorted.bam ./alignment/unsorted_bam/Nfxl1_GM12878_input_rep1_unsorted.bam
sambamba-0.8.0 view -h -t 4 -f bam -F "[XS] == null and proper_pair and not (unmapped or mate_is_unmapped) and not duplicate" ./alignment/sorted_bam/Nfxl1_GM12878_input_rep1_sorted.bam > ./alignment/filtered_bam/Nfxl1_GM12878_input_rep1_aln.bam
sambamba-0.8.0 sort -t 4 -o ./alignment/sorted_filtered_bam/Nfxl1_GM12878_input_rep1_aln_sorted.bam ./alignment/filtered_bam/Nfxl1_GM12878_input_rep1_aln.bam

```

There isn’t much difference in the commands for single end and paired end reads except for the final filtering step. ‘XS’ is a tag used by Bowtie2 that gives an alignment score for the second-best alignment, and it is only present if the read is aligned more than once.

Also the ‘t’ parameter refers to the number of threads available for computing and you can adjust that to suit your system.

You can run the same command for each sample or use the script samples_postalign.sh from **the scripts folder in the repository.**
 
--------------------------------------------------------------------

### 5. Peak Calling of Aligned ChIP-Seq Reads <a name="peak_calling"></a>

Peak calling is the most important step of ChIP-Seq analysis that helps identify genomic regions to which transcription factors or histone marks are bound. Identifying such genomic regions helps us understand which genes are being regulated.

The binding sites to which the proteins are bound appear as areas with high read density and are known as peaks. Sequencing of DNA from ChIP samples is performed randomly from either end and it does not cover the entire length of the enriched DNA fragment.

The 5’ end of the selected fragments from groups on the positive and negative strand centered around the binding site. (See the figure below). These binding patterns form a characteristic bimodal distribution or paired peaks.

<p> </br> </p>
<p align="center">
<img src="https://github.com/ShrutiBaikerikar/ChipSeq_Tutorial/blob/main/images/chip_seq_bimodal_peaks_image_1.png" width="400" height=500 alt="Bimodal Binding Pattern of ChIP-Seq Reads image"/>
</p>

<p align="center">
     <b> Characteristic BiModal Binding Pattern of ChIP-Seq Reads </b>
</p> 

*Citation: Wilbanks EG, Facciotti MT (2010) Evaluation of Algorithm Performance in ChIP-Seq Peak Detection. PLoS ONE 5(7): e11471. https://doi.org/10.1371/journal.pone.0011471*

This characteristic shape obtained from ChIP-Seq samples is assessed using statistical measures and compared against background (input samples) to confirm whether the site of enrichment is really a potential binding site.

<p> </br> </p>
<p align="center">
<img src="https://github.com/ShrutiBaikerikar/ChipSeq_Tutorial/blob/main/images/chip_seq_bimodal_peaks_image_2.png" width="600" height=800 alt="Different Signal Types from Chip-Seq Data image"/>
</p>

<p align="center">
     <b> Types of Peaks or Signals from different ChIP-Seq data </b>
</p> 

*Citation: Pepke S, Wold B, Mortazavi A. Computation for ChIP-seq and RNA-seq studies. Nat Methods. 2009;6(11 Suppl):S22-S32. doi:10.1038/nmeth.1371*

While ChIP-seq samples have a characteristic binding pattern, the type of signals differ. In the above image: a] Narrow peak for transcription factor CTCF, b] Mixed signal for RNA Pol III c] Broad peaks for histone mark H3K27me3

Here is a primer as to why we observe different signal types:
- Narrow Peaks/Sharp signal: Transcription factors generally bind to very specific DNA sequence motifs. Thus the enriched fragments are centered around the motif leading to sharp, peaky regions.
- Broad Peaks/Broad signal: Histone marks are spread over multiple nucleosomes that are loosely positioned on the DNA sequence; in other words, histone marks cover large gene bodies. Hence the signal appears as broad regions of enrichment and are several kilobases in size.
- Mixed signal: Mixed binding signals are mostly observed for RNA polymerase i.e. PolII/III show broad and sharp peaks. Narrow peaks generally suggest binding at the promotor while broad peaks correspond to transcription elongation.

#### 5a. Peak Calling with MACS2 - Theory

Primarily, the characteristic bimodal distribution of ChIP-seq samples assists in peak calling. Peak callers generally follow this basic framework:
- slide a window along the genome to identify regions of enrichment
- calculate enrichment of reads in ChIP samples vs input samples
- calculate a significance score corrected for multiple testing.

In this tutorial, we will be using the peak caller MACS2 or Model-Based Analysis of ChIP-Seq. MACS2 is written in Python and you can download it from https://github.com/macs3-project/MACS

Download MACS2, install it, and add it to your path. MACS2 requires Python to function so ensure you have Python installed and added to your path as well.

MACS2 performs peak calling by taking advantage of the bimodal distribution of ChIP-Seq reads. It takes into account information about sequencing read position, its orientation as well as the influence of genome complexity on the signals.

MACS2 works for both single-end and paired-end reads. It works for ChIP sample alone or along with an input sample that increases the specificity of peak calls.

Let’s go over the steps followed MACS2 algorithm in detail.

**Step 1: Removing duplicates**

A good sequencing depth can help identify all true binding sites in the sample. But sequencing is susceptible to errors and these can give rise to duplicates.

Here reads with the same start position are considered duplicates. There can be two types of duplicates:
- Bad duplicates: These can arise from PCR overamplification or biases in sequence library preparation. Also blacklisted or repeat regions in the genome with excessively high signals can contribute to duplicates. Such duplicates add noise to final peak calls and hence need to be eliminated.
- Good duplicates: Increasing sequencing depth can result in valid biological duplicates or reads that map at the same position and are biologically relevant. Removal of biological duplicates can lead to under-estimated ChIP-Seq signals.

MACS2 deals with duplicates in several ways via –KEEP-DUP command. The default auto option makes MACS calculate the maximum tags at the exact same location based on binomial distribution using 1e-5 as the p-value cutoff.
The all option keeps every tag.

You can specify an integer and accordingly MACS will keep at the most that many tags at the same location. By default, it keeps only one read at the same location. 

**Step 2: Modeling the shift size**

As mentioned previously, ChIP sequencing is performed randomly from either end. This results in reads that represent the ends of the fragments instead of the actual protein-DNA binding sites.

The 5’end of the reads are centered around the binding site both on the positive and negative strand thus resulting in a bimodal distribution. For single-end reads, the distance between the forward and reverse strand distribution is used to estimate the fragment size.

In the case of paired-end reads, MACS2 estimates fragment size from the insert size of the paired reads.

Estimation of fragment size helps shift the 5’end fragments of Chip-Seq samples towards the 3’ direction to better estimate the true binding site. The distance by which the 5’end ends are to be moved is referred to as shift size.

To model the shift size, MACS2 scans the genome with a window of length 2 x bandwidth and identifies regions with reads more than the mfold enriched relative to random tag distribution. 
Bandwidth refers to the sonication size or the size of DNA fragments that are sheared after immunoprecipitation. Mfold is high-confidence fold-enrichment.

Then MACS2 randomly samples 1000 of these high-quality peaks, separates the positive and negative strands, and aligns them at the midpoint between their centers. (See the figure below)

<p> </br> </p>
<p align="center">
<img src="https://github.com/ShrutiBaikerikar/ChipSeq_Tutorial/blob/main/images/chip_seq_bimodal_peaks_image_3.png" width="400" height=500 alt="Peak Calling with MACS2 image"/>
</p>

<p align="center">
     <b> Theoretical Visualization of MACS2's modeling of shift size  </b>
</p> 

*Citation: Bardet AF, Steinmann J, Bafna S, Knoblich JA, Zeitlinger J, Stark A. Identification of transcription factor binding sites from ChIP-seq data at high resolution. Bioinformatics. 2013;29(21):2705-2713. doi:10.1093/bioinformatics/btt470*

The distance between the modes of the forward and reverse strand peaks in the alignment is defined as 'd', and MACS2 shifts all the tags by d/2 toward the 3' ends to the most likely protein-DNA interaction sites.

**Step 3: Peak detection**

Before peak detection, for experiments with control samples, MACS2 linearly scales the total control read count to be the same as the total ChIP read count.

The read distribution of ChIP-Seq data is modelled as the Poisson distribution. 

Technically, the Poisson distribution helps predict the probability of a given number of events happening in a fixed interval of time. But here we use Poisson distribution to model read distribution in ChIP-Seq because:
- We are dealing with count data.
- The distribution is applicable for data with zero or positive integers.
- And the variance is a function of mean.

For Poisson distribution, variance is equal to mean and is sometimes written as lambda. It refers to the rate of the event.

However, count data such as read distribution tends to have variance larger than the mean which is referred to as overdispersion. This is generally modelled by negative binomial distribution and there are a couple of peak callers that implement the same.

MACS2 uses Poisson distribution to model tag distribution along the genome but with a slight modification.

First, it shifts every tag by d/2, then it slides a 2d window across the genome to find candidate peaks with a significant tag enrichment. The significance is measured by p-value from Poisson distribution with λBG parameter.

λBG (BG refers to background genome) captures both the mean and variance of the distribution. It is the expected number of reads in that window.

To calculate λBG , MACS2 requires the effective genome size. Effective genome size is defined as the genome size which can be sequenced. Due to the repetitive features on the chromosomes, the actual mappable genome size is smaller than the original size, about 90% or 70% of the genome size.

For the human genome, it is 2.7e9. MACS2 provides pre-computed values for human, mouse, worm and fly genomes.

Now coming to the modification that MACS2 uses to model read distribution from ChIP-Seq. In control samples, read distributions have fluctuations and biases arising from DNA amplification and sequencing, genome copy number variation ,and chromatin structure.

Therefore instead of using a uniform λBG estimated from the whole genome, MACS2 uses a dynamic parameter, λlocal, defined for each candidate peak as:

λlocal = max(λBG, [λ1k,] λ5k, λ10k)
where λ1k, λ5k ,and λ10k are λ estimated from the 1 kb, 5 kb ,or 10 kb window centered at the peak location in the control sample, or the ChIP-Seq sample when a control sample is not available (in which case λ1k is not used)

λlocal captures local biases. This makes MACS2 robust against the occasional low read count in certain relevant genomic regions. It also removes potential false positives due to local biases (that is, peaks significantly under λBG, but not under λlocal).
MACS2 calculates the p-value for every candidate peak using λlocal. 

Overlapping enriched peaks are merged, and each read position is extended d bases from its center. The location with the highest fragment pileup, referred to as the summit, is predicted as the precise binding location.

Candidate peaks with p-values below a user-defined threshold p-value (default 10-5) are called, and the ratio between the ChIP-Seq tag count and λlocal is reported as the fold_enrichment.

**Step 4: Multiple testing correction**

Each peak is considered as an independent test and since we are detecting thousands of significant peaks, we run into multiple testing problem.

If you run a single hypothesis test, there’s a small chance (alpha = 5%) that you’ll get a false significant result. If you run thousands of tests and maintain the same alpha value, then the number of false alarms increases dramatically.
This can be avoided with multiple testing correction.

Initially, MACS estimated the False Discovery Rate by finding ChIP peaks over control and control peaks over ChIP (sample swap). FDR was defined as Number of control peaks / Number of ChIP peaks.
Now MACS2 corrects the p-values with Benjamini-Hochberg correction.

#### 5b. Running MACS2 for calling Broad and Narrow Peaks

We have alignment files for histone marks and transcription factors, so we will run separate commands to call broadPeaks and narrowPeaks with MACS2.

Also, we have single-end reads (for histones) and paired end reads (for TFs) so we will specify separate file formats BAM and BAMPE respectively while calling peaks with MACS.

To learn about all the possible commands for executing MACS2, please refer to https://pypi.org/project/MACS2/

If you have installed and added MACS2 to your path, run the following commands:

**For Histone Marks**

```bash

macs2 callpeak -t ./alignment/filtered_bam/H3K4me1_GM12878_rep1_aln.bam -c ./alignment/filtered_bam/H3K4me1_GM12878_input_rep1_aln.bam -f BAM --broad -g hs --broad-cutoff 0.1 --keep-dup auto --outdir ./macs2/output -n H3K4me1_GM12878_rep1

```

Details of the above-mentioned command are:
- callpeak : Main MACS2 Function to call peaks from alignment results.
- t : Immunoprecipitation File. In this case, it is H3K4me1_GM12878_rep1_aln.bam
- c: Control/Input File. In this case, it is H3K4me1_GM12878_input_rep1_aln.bam
- f : Format of the input file. In this case, it is single-end BAM file, hence f is specified as BAM
- broad: To identify broad peaks, MACS tries to combine broad regions by putting nearby highly enriched regions into a broad region with loose cutoff.
- broad-cutoff: Cutoff for broad region.
- g: The mappable genome size or effective genome size. In this case, the default hs -- 2.7e9 is applicable for human genome.
- keep-dup: It controls the MACS behavior towards the duplicate tags at the exact same location. The default auto option makes MACS calculate the maximum tags at the exact same location based on binomial distribution using 1e-5 as p-value cutoff.
- outdir: Specify the directory to store macs2 output files
- n: Specify the name string of the experiment that MACS can use while creating output files.

**For Transcription Factors**

```bash

macs2 callpeak -t ./alignment/filtered_bam/Nfxl1_GM12878_rep1_aln.bam -c ./alignment/filtered_bam/Nfxl1_GM12878_input_rep1_aln.bam -f BAMPE -g hs --keep-dup auto --outdir ./macs2/output -n Nfxl1_GM12878_rep1 -B

```	 
Details of the above-mentioned command are:
- callpeak : Main MACS2 Function to call peaks from alignment results.
- t : Immunoprecipitation File. In this case, it is H3K4me1_GM12878_rep1_aln.bam
- c: Control/Input File. In this case, it is H3K4me1_GM12878_input_rep1_aln.bam
- f : Format of the input file. In this case, it is paired-end BAM file, hence f is specified as BAMPE
- g: The mappable genome size or effective genome size. In this case, the default hs -- 2.7e9 is applicable for human genome.
- keep-dup: It controls the MACS behavior towards the duplicate tags at the exact same location. The default auto option makes MACS calculate the maximum tags at the exact same location based on binomial distribution using 1e-5 as p-value cutoff.
- outdir: Specify the directory to store macs2 output files
- n: Specify the name string of the experiment that MACS can use while creating output files.
- B: With this flag specified, MACS stores the fragment pileup, control lambda in bedGraph files.

You can run the same command for each sample and its replicate or use the script samples_peakcalling.sh from **the scripts folder in the repository.**

#### 5c. Understanding MACS2 output files

Before, we get to the output files obtained from MACS2, lets understand the different file formats used in ChIP-Seq.

**BED Format**
It is a tab-delimited file that contains data about genomic coordinates of various genomic “features” such as exon, UTRs, etc. It has 0-based coordinates.
The first three fields/columns in each feature line are required:
- chr: chromosome name/ID
- start: start position of the feature
- end: end position of the feature
There are nine additional fields that are optional.

**BEDGraph Format**
It is used to display continuous-valued data in a track format; It helps display density or coverage information. It stores data in the original format.
It is based on the BED format but the ‘score’ is stored in column 4 instead of 5 and track lines are included.

**Wiggle Format**
It is similar to BEDGraph but the data is compressed, data elements are stored in bins of equal size and it has 1-based coordinates. 
It associates a floating point number with positions in the genome, which is plotted on the track’s vertical axis to create a wiggly line.

**bigWig Format**
It is an indexed binary format derived from Wiggle format.

Now coming to the output files generated by MACS:
1. NAME_peaks.xls is a tabular file which contains information about called peaks. 
2. NAME_peaks.narrowPeak is BED6+4 format file which contains the peak locations together with peak summit, p-value, and q-value. 
3. NAME_summits.bed is in BED format, which contains the peak summits locations for every peak. The 5th column in this file is the same as what is in the narrowPeak file. 
4. NAME_peaks.broadPeak is in BED6+3 format which is similar to the narrowPeak file, except for missing the 10th column for annotating peak summits. This file and the gappedPeak file will only be available when --broad is enabled. 
5. NAME_peaks.gappedPeak is in BED12+3 format which contains both the broad region and narrow peaks. 
6. NAME_model.r is an R script which you can use to produce a PDF image of the model based on your data. This file is not available for paired-end reads or BAMPE format as fragment length is determined by the insert size of paired end reads and not from the bimodal distribution of plus and minus strand reads.
7. NAME_treat_pileup.bdg and NAME_control_lambda.bdg files are in bedGraph format The NAME_treat_pielup.bdg contains the pileup signals from ChIP/treatment sample. The NAME_control_lambda.bdg contains local biases estimated for each genomic location from the control sample, or from treatment sample when the control sample is absent. 

As per the commands that we have executed above, we should get 4 files each for histone marks( RScript, gappedPeak, broadPeak and .xls file) and 5 files each for transcription factor (narrowPeak, .xls file, summits.bed file, control_lambda.bdg, treat_pileup.bdg).

 
------------------------------------------------------------------------

### 6. ChIP-specific Quality Control <a name="chip_qc"></a>

ChIP-seq data quality is again examined via specific metrics after genomic alignment. ChIPQC is a R/Bioconductor package that provides a comprehensive ChIP-seq-specific quality control analysis.

When conducting ChIP-QC analysis, it is advisable to use RStudio. Please install RStudio from https://www.rstudio.com/

To conduct ChIP-QC analysis, you will also need the blacklist region file for the human genome (GRCh38).  Please download the file ENCFF356LFX from https://www.encodeproject.org/annotations/ENCSR636HFF/

Once downloaded, rename and save the file as GRCh38_unified_blacklist.bed in the blacklist folder in the data folder.

Blacklist regions comprise genomic regions that tend to have a very high ratio of multi-mapping to unique mapping reads and high variance in mappability. They include repeat elements such as satellite, centromeric and telomeric repeats.

Such regions have anomalous or excessively high signals in next-generation sequencing experiments independent of cell lines or experiments. The removal of the ENCODE blacklist is an essential quality measure when analyzing functional genomics data.

Before we run the ChIP-specific analysis, we need to create a sample sheet contains the metadata information for our dataset. It includes the following columns:
* SampleID: Identifier string for sample
* Tissue, Factor, Condition: Identifier strings for up to three different experimental parameters such as tissue/cell line, factor indicates the name of the protein you are studying (histone mark or TF) and condition would be treatment/no treatment ,etc. In our case, we simply want to differentiate between our three sample sets (H3K4me1_GM12878, Nfxl1_GM12878, Nfxl1_K562) and that’s what we have set our Condition column as.
If you don’t have information, then set values to NA.
* Replicate: Replicate number of each sample
* bamReads: file path for BAM file containing aligned reads for ChIP sample
* ControlID: an identifier string for the control sample
* bamControl: file path for bam file containing aligned reads for the control sample
* Peaks: path for file containing peaks for the sample
* PeakCaller: Identifier string for peak caller used. In our case, we have set it to “bed”.

Both the meta_samples_chipqc.csv file and the chipqc.R script has been provided in the **tutorial_data and scripts folder in the repository respectively**. Please change the file path in the meta_samples_chipqc.csv to suit your system.

For this analysis, we will use the chipseq_r folder which contains subdirectories bam and peaks.
Before we begin, copy the (_aln_sorted.bam) BAM files and their respective indices from ./align/sorted_filtered_bam to ./chipseq_r/bam folder. Also copy the .narrowPeak and broadPeak files for all samples from ./macs2/output to ./chipseq_r/peaks folder.

Start a session in RStudio and set chipseq_r as your working directory. Create a new R Script ‘chipqc.R’ and run the following commands:

**COMMANDS**

```r

BiocManager::install(“BiocParallel”)
BiocManager::install(“ChIPQC”)
BiocManager::install(“TxDb.Hsapiens.UCSC.hg38.knownGene”)

library(BiocParallel)
library(ChIPQC)
library(TxDb.Hsapiens.UCSC.hg38.knownGene)

#load meta data of the samples
#Change the path to the file meta_samples_chipqc.csv to suit your system
samples <- read.csv('chip_seq_analysis/chipseq_r/meta_samples_chipqc.csv')
samples

#If you are using a Windows PC, please run this command
register(SerialParam())

## Create ChIPQC object [This step takes some time.]
chipObj <- ChIPQC(samples, annotation="hg38", consensus = TRUE, blacklist = "chip_seq_analysis/data/blacklist/GRCh38_unified_blacklist.bed" ) 

## Create ChIPQC report
ChIPQCreport(chipObj, reportName="ChIP QC report", reportFolder="ChIPQCreport", facetBy = "Condition", lineBy = "Replicate")

```
Running these commands will generate a ChipQC Folder containing around 9 plots and HTML file containing the ChIPQC results. Lets examine these further.

#### 6a. Understanding ChIPQC Report

<p> </br> </p>
<p align="center">
<img src="https://github.com/ShrutiBaikerikar/ChipSeq_Tutorial/blob/main/images/chipqc_image_1.png" width="800" height=400 alt="Summary of ChIPQC of all samples image"/>
</p>

<p align="center">
     <b>ChipQC Report: Short summary of QC checks for all samples </b>
</p> 

When you open the html file that contains the ChIPQC report, you first get a short summary about the quality control checks for your samples.

In this summary, various metrics are covered:
- Reads - Number of sample reads within analysed chromosomes.
- Dup% - Percentage of MapQ filter passing reads marked as duplicates
- FragLen - Estimated fragment length by cross-coverage method
- SSD - SSD score 
- FragLenCC - Cross-Coverage score at the fragment length
- RelativeCC - Cross-coverage score at the fragment length over Cross-coverage at the read length
- RIP% - Percentage of reads within the peaks
- RIBL% - Percentage of reads within the Blacklist regions

These metrics helps us understand the distribution of the signal in the enriched regions, in the potential artefact regions (blacklist regions) as well as across the whole genome.

**SSD Score**

<p> </br> </p>
<p align="center">
<img src="https://github.com/ShrutiBaikerikar/ChipSeq_Tutorial/blob/main/images/chipqc_image_2.png" width="800" height=400 alt="ChIPQC: SSD of all samples image"/>
</p>

<p align="center">
     <b>ChipQC Report: Standardised Standard Deviation of the coverage plot for all samples </b>
</p> 

SSD Score stands for Standardised Standard Deviation of the coverage. It the standard deviation of coverage normalised to the total number of reads.

A higher proportion of genomic positions at greater depths is indicative of successful ChIP-Seq samples. A high SSD score means that there are regions with high signal in the sample but this can be strongly affected by artefacts.

In our dataset, Nfxl1 samples have a higher SSD score than H3K4me1. The coverage plot also indicates better pileup at higher sequencing depth or heavier tail which is indicative of good enrichment.

**RiP%**

<p> </br> </p>
<p align="center">
<img src="https://github.com/ShrutiBaikerikar/ChipSeq_Tutorial/blob/main/images/chipqc_image_3.png" width="800" height=400 alt="ChIPQC: Percentage of Reads in Peaks for all samples image"/>
</p>

<p align="center">
     <b>ChipQC Report: Percentage of Reads in Peaks for all samples </b>
</p> 


RiP% (also known as FRiP%) is the percentage of reads in peaks. This indicates how enriched the sample is or how successful the immunoprecipitation reaction is.

It is generally suggested that in a good ChIP-seq library, FRIP should be >5% for a TF and >25-30% for broad histone modifications or RNA polymerase II. But various other sources comment that there are known examples of good datasets with FRiP < 1%.

In our samples, RiP% is greater for H3K4me1 samples compared to Nfxl1 which is obvious since histone marks have broad peaks that span over larger genomic regions.

Also the plot shows that ChIP samples are more enriched than input samples.

**RiBL%**

<p> </br> </p>
<p align="center">
<img src="https://github.com/ShrutiBaikerikar/ChipSeq_Tutorial/blob/main/images/chipqc_image_4.png" width="800" height=400 alt="ChIPQC: Percentage of Reads in Blacklist for all samples image"/>
</p>

<p align="center">
     <b>ChipQC Report: ChIPQC: Percentage of Reads in Blacklist for all samples </b>
</p> 


RiBL% is the percentage of reads in the blacklisted regions. Blacklisted regions are genomic regions mostly comprising of repeats that have excessively high signals in NGS experiments and add noise to the true signal.
It is best to have low RiBL percentages. For our dataset, we observe the same in the report summary as well as in the plot.

**Enrichment in expected genomic regions**

<p> </br> </p>
<p align="center">
<img src="https://github.com/ShrutiBaikerikar/ChipSeq_Tutorial/blob/main/images/chipqc_image_5.png" width="800" height=400 alt="ChIPQC: Enrichment in expected genomic regions for all samples image"/>
</p>

<p align="center">
     <b>ChipQC Report: Plot of Enrichment in Expected genomic regions for all samples </b>
</p> 

Most factors will be enriched in one or a few of the annotated genomic regions: promoters, exons, introns, 5'or 3' untranslated regions (UTRs). 

The relative enrichment of ChIP-seq reads within these regions compared to the features’ genomic proportion indicates whether the expected genomic distribution was captured.
In our dataset, we see an enrichment in ‘Promoters’ and ‘All 5’utrs’ as well as a mild enrichment in ‘All cds’ for both histone marks and TFs.

**FragL and RelCC**

<p> </br> </p>
<p align="center">
<img src="https://github.com/ShrutiBaikerikar/ChipSeq_Tutorial/blob/main/images/chipqc_image_6.png" width="800" height=400 alt="ChIPQC: Cross Coverage Plot for all samples image"/>
</p>

<p align="center">
     <b>ChipQC Report: Cross Coverage Plot for all samples </b>
</p> 

FragLength (Fragment Length Cross Coverage) and RelCC (Relative Strand Cross Correlation Coefficient) are statistics related to the strand cross-correlation.

As mentioned previously, some ChIP-Seq reads map to the forward strand forming a normal distribution upstream of the TF binding site, while those mapping to the reverse strand mirror the same distribution downstream. 
In other words, we have a bimodal distribution around the binding site. The distance between the modes of these distributions equals the average length of the enriched DNA fragments.

RelCC values larger than 1 for all ChIP samples suggest good enrichment. We observe this for our dataset.

The FragL values should be roughly the same as the fragment length picked during the size selection step in library preparation. This generally ranges from 100-600bp.

In a cross-correlation analysis, the reverse reads are shifted towards the forward reads. The cross-correlation metric is computed as the Pearson’s linear correlation between coverage for each base on the minus strand and the plus strands, by systematically shifting minus strand by k base pairs at a time. 
This shift is performed over and over again to obtain the correlations for a given area of the genome The cross-correlation between forward and reverse reads should be highest at the shift which corresponds to the average fragment size.

ChIPQC calculates a slightly different metric called Cross-Coverage Score. It is based on the idea that shifting the reverse reads forward should reduce the total number of base pairs covered in the genome. 
The lowest coverage should be obtained when the reverse reads are shifted by the fragment length, and thus cross-coverage is highest at this point.

Thus, in the CC Plot the highest peak indicates Fragment Length while the second peak indicates read length or results from blacklisted or overlapped regions. The Fragment Length peak should be higher than read length peak.

We have got high Fragment Length peaks for our samples in the plot (especially in comparison to input samples) and good RelCC scores. But the read length peak is not obvious. The peaks are not clear for histone marks as they are broad peaks covering longer genomic areas. 

<p> </br> </p>
<p align="center">
<img src="https://github.com/ShrutiBaikerikar/ChipSeq_Tutorial/blob/main/images/chipqc_image_7.png" width="800" height=400 alt="Example of Cross Coverage Plots image"/>
</p>

<p align="center">
     <b>Example of Cross Coverage Plots for different ChIP-Seq experiments </b>
</p> 

*Citation: Landt SG, Marinov GK, Kundaje A, et al. ChIP-seq guidelines and practices of the ENCODE and modENCODE consortia. Genome Res. 2012;22(9):1813-1831. doi:10.1101/gr.136184.111*

This image shows a series of CC Plot indicative of successful, marginal or failed ChIP-Seq experiments and gives us a general idea of what our CC plots should look like.

**Peak Profile**

<p> </br> </p>
<p align="center">
<img src="https://github.com/ShrutiBaikerikar/ChipSeq_Tutorial/blob/main/images/chipqc_image_8.png" width="800" height=400 alt="ChIPQC: Peak profile for all samples image"/>
</p>

<p align="center">
     <b>ChipQC Report: Peak profile for all samples </b>
</p> 

This plot represents the mean read depth across and around peaks. 
The shape of these profiles can vary depending on what type of mark is being studied – transcription factor, histone mark, or polymerase; but similar distinct marks differing from input profile indicate successful ChIPs.

Now that we have assessed our samples with ChIP specific quality control measures and also evaluated the presence of blacklisted regions, we can filter the peaks corresponding to the blacklisted regions.

This step can be done before peak calling and is usually done immediately after peak calling before ChIPQC, but to demonstrate the percentage of blacklisted regions in our samples we chose to filter the peaklist post ChIPQC.

We are using BedTools to filter out blacklisted regions. You can download this from https://bedtools.readthedocs.io/en/latest/index.html. After installing and adding it to your path, run the following commands.

```bash

bedtools intersect -v -a ./macs2/output/H3K4me1_GM12878_rep1_peaks.broadPeak -b ./data/blacklist/GRCh38_unified_blacklist.bed > ./macs2/output/H3K4me1_GM12878_rep1_blfilt_peaks.broadPeak

bedtools intersect -v -a ./macs2/output/Nfxl1_GM12878_rep1_peaks.narrowPeak -b ./data/blacklist/GRCh38_unified_blacklist.bed > ./macs2/output/Nfxl1_GM12878_rep1_blfilt_peaks.narrowPeak

```

Here, I have shown only one example of .narrowPeak and .broadPeak file. Remember to run these commands for each of the .narrowPeak and .broadPeak files obtained after running MACS2.
These blacklist filtered files are available in **the tutorial_data folder in the repository.**

--------------------------------------------------------------------

### 7. Handling Replicates with IDR and BedTools <a name="replicates"></a>

Comparative analysis of ChIP-Sequencing data helps compare the binding of the protein of interest either between replicates or under different physiological conditions. 

Every experimental design requires at least 2-3 biological replicates to overcome the amount of variability in the data that one would get with only one sample. With more than one replicate you could observe a common trend in your data more frequently for the given condition thus strengthening your findings.

Technically, multiple replicates measuring the same biological conditions should have high consistency. But that is not always the case. One reason for this is that every comparison in ChIP-Seq data analysis is strongly influenced by the thresholds that were applied during the peak calling step.

To evaluate consistency between replicates in ChIP-Seq, several methods have been designed such as finding overlapping peaks and assessing the difference in binding regions and complex statistical methods that evaluate reproducibility between replicates.

#### 7a. Handling Replicates with IDR

IDR stands for Irreproducibility Discovery Rate. The ENCODE project developed the irreproducible discovery rate (IDR), as a metric to identify genuine peaks based on their reproducibility across replicates.

IDR begins with generating a peak list with a lenient threshold that will contain both genuine peaks and noise. The peaks in the list are ranked, e.g. on their fold enrichment or significance, these ranks will correlate well for genuine peaks, while the noise will show no correlation. 

Further, IDR uses a statistical method to find the point in the curve at which the heterogeneity of the association between replicate ranks sharply increases. This decay point, from which no consistency in correspondence is observed, serves as an internal indicator of the change from signal to noise.

As a metric to evaluate consistency between replicates, IDR has several advantages:
- It does not take into account initial cut-offs used during peak calling as they tend to be different for different peak callers. This means that you can compare peak lists obtained from different peak callers.  
- It does not depend on specific thresholds, so all peaks are considered.

IDR is not recommended for broad signals from histone marks as there is an ambiguous overlap of broad peaks.

The ENCODE3 guidelines suggest that IDR should be conducted to compare the following pairs:
- True replicates
- Pooled pseudo-replicates
- Self-pseudo-replicates for each replicate

The results of these comparisons are further evaluated to obtain the final peak call set for downstream analysis. In this tutorial, we will be comparing only True replicates.

**Running IDR**

You can download IDR from https://github.com/nboley/idr This one is written in Python and is run via the Bash shell. We will be using this version for the tutorial.

Additionally IDR is written in R as well and you can learn more about it here https://www.encodeproject.org/software/idr/

After installing IDR and adding to your path, run the following commands:

```bash

#Call peaks with relaxed thresholds on the first replicate
macs2 callpeak -t ./alignment/filtered_bam/Nfxl1_GM12878_rep1_aln.bam -c ./alignment/filtered_bam/Nfxl1_GM12878_input_rep1_aln.bam -f BAMPE -g hs -p 1e-3 --keep-dup auto --outdir ./macs2/output2 -n Nfxl1_GM12878_rep1_relaxed -B

#Count peaks called with relaxed threshold in replicate 1
wc -l ./macs2/output2/Nfxl1_GM12878_rep1_relaxed_peaks.narrowPeak

#Call peaks with relaxed threshold on the second replicate
macs2 callpeak -t ./alignment/filtered_bam/Nfxl1_GM12878_rep2_aln.bam -c ./alignment/filtered_bam/Nfxl1_GM12878_input_rep2_aln.bam -f BAMPE -g hs -p 1e-3 --keep-dup auto --outdir ./macs2/output2 -n Nfxl1_GM12878_rep2_relaxed -B

#Count the peaks called with relaxed threshold in the second replicate
wc -l ./macs2/output2/Nfxl1_GM12878_rep2_relaxed_peaks.narrowPeak

#Create pooled ChIP sample from replicates
samtools merge --threads 4 ./alignment/pooled_bam/Nfxl1_GM12878_pooled.bam ./alignment/filtered_bam/Nfxl1_GM12878_rep1_aln.bam ./alignment/filtered_bam/Nfxl1_GM12878_rep2_aln.bam 

#Create pooled input sample from inputs for each replicate
samtools merge --threads 4 ./alignment/pooled_bam/Nfxl1_GM12878_input_pooled.bam  ./alignment/filtered_bam/Nfxl1_GM12878_input_rep1_aln.bam ./alignment/filtered_bam/Nfxl1_GM12878_input_rep2_aln.bam 

#Call peaks with relaxed threshold on pooled sample
macs2 callpeak -t ./alignment/pooled_bam/Nfxl1_GM12878_pooled.bam -c ./alignment/pooled_bam/Nfxl1_GM12878_input_pooled.bam -f BAMPE -g hs -p 1e-3 --keep-dup auto --outdir ./macs2/output2 -n Nfxl1_GM12878_pooled_relaxed -B

#Count peaks called with relaxed threshold on pooled sample
wc -l ./macs2/output2/Nfxl1_GM12878_pooled_relaxed_peaks.narrowPeak

```
The above code focuses only on samples of transcription factor Nfxl1 in cell line GM12878. You have to run the same commands for Nfxl1 in cell line K562. Also remember to change the 'threads' parameter to suit your system.

Now, that we have all our peak list files (with relaxed threshold) ready, we can conduct an IDR analysis on the replicates. 

```bash

#Run IDR analysis on replicates with relaxed threshold
idr --samples ./macs2/output2/Nfxl1_GM12878_rep1_relaxed_peaks.narrowPeak ./macs2/output2/Nfxl1_GM12878_rep2_relaxed_peaks.narrowPeak --peak-list ./macs2/output2/Nfxl1_GM12878_pooled_relaxed_peaks.narrowPeak --input-file-type narrowPeak --rank p.value --output-file ./idr/Nfxl1_GM12878_idr_output --soft-idr-threshold 0.05 --plot --use-best-multisummit-IDR

```
Run the same command for the files generated from peak calling with relaxed thresholds on Nfxl1 in K562 cell line.

To run the IDR analysis, we have used the following parameters:
- samples: list of narrowPeak files for each replicate that have been called with relaxed thresholds
- peak-list: narrowPeak file containing peaks called with relaxed threshold on pooled samples. For every peak in this file a single peak from each replicate is chosen that overlaps this peak.
- input-file-type: in our case, narrowPeak as we are examining transcription factors
- rank: the statistical measure used to rank peaks; it could be signal.value, p.value, q.value. 
- output-file: path and name of output file
- soft-idr-threshold: Report statistics for peaks with a global idr below this value but return all peaks with an idr below --idr. --idr is IDR_THRESHOLD which only return peaks with a global idr threshold below this value
- plot: plots the results

Now that we have conducted IDR analysis, lets extract our consensus set of peaks for further downstream analysis. We will use only those peaks that cross the IDR threshold of 0.05.

```bash

#Defining the IDR threshold
IDR_THRESH_TRANSFORMED=$(awk -v p=0.05 'BEGIN{print -log(p)/log(10)}')

#Getting peaks that pass the IDR threshold
awk 'BEGIN{OFS="\t"} $12>='"${IDR_THRESH_TRANSFORMED}"' {print $1,$2,$3,$4,$5,$6,$7,$8,$9,$10}' ./idr/Nfxl1_GM12878_idr_output | sort | uniq | sort -k7n,7n > ./idr/Nfxl1_GM12878_idr0.05.narrowPeak

#Counting the peaks that crossed the IDR threshold
wc -l ./idr/Nfxl1_GM12878_idr0.05.narrowPeak

#Filter using a blacklist and remove non-standard chromosome names
bedtools intersect -v -a ./idr/Nfxl1_GM12878_idr0.05.narrowPeak -b ./data/blacklist/GRCh38_unified_blacklist.bed | grep -P 'chr[\dXY]+[ \t]' > ./idr/Nfxl1_GM12878_idr0.05_filt.narrowPeak

#Counting the peaks in the blacklist filtered file
wc -l ./idr/Nfxl1_GM12878_idr0.05_filt.narrowPeak 

```

**Output Files of IDR**
The output file mimics the input file type which in our case is the .narrowPeak file. Apart from that, IDR also produces plots pertaining to the samples.

Some interesting features about the narrowPeak file obtained from IDR include:
- Column no 5: Contains the scaled IDR value, min(int(log2(-125IDR), 1000). e.g. peaks with an IDR of 0 have a score of 1000, idr 0.05 have a score of int(-125log2(0.05)) = 540, and idr 1.0 has a score of 0.
- Column no 11: Contains Local IDR [-log10(Local IDR value)]. It is similar to the posterior probability of a peak belonging to the irreproducible noise component.
- Column no 12: Contains Global IDR [-log10(Global IDR value)]. The global IDR is the value used to calculate the scaled IDR number in column 5; it is analogous to a multiple hypothesis correction on a p-value to compute an FDR.

<p> </br> </p>
<p align="center">
<img src="https://github.com/ShrutiBaikerikar/ChipSeq_Tutorial/blob/main/images/idr_plots_image1.png" width="800" height=400 alt="IDR output Plots for all samples image"/>
</p>

<p align="center">
     <b>IDR Output plots for both samples Nfxl1 in GM12878 and Nfxl1 in K562 cell line </b>

These are the plots obtained for our transcription factor samples in the dataset: Nfxl1_GM12878 and Nfxl1_K562.

Here are the details of the plots as mentioned by the creators of IDR:
- Upper Left: Replicate 1 peak ranks versus replicate 2 peak ranks - peaks that do not pass the specified IDR threshold are coloured red.
- Upper Right: Replicate 1 log10 peak scores versus replicate 2 log10 peak scores - peaks that do not pass the specified IDR threshold are coloured red.
- Bottom Row: Peaks rank versus IDR scores are plotted in black for each replicate. The overlayed boxplots display the distribution of IDR values in each 5% quantile. 


#### 7b. Overlapping of Peaks with BedTools

Overlap of peak regions is a straightforward approach to compare peaks across the samples and develop a consensus peak set. It helps us understand whether a peak is found in both samples and gives us a rough idea of the similarity in peak regions between different samples.

However, this approach underestimates the similarity because it does not take into account the different thresholds used during peak calling.

According to the ENCODE3 guidelines, IDR analysis is not recommended for histone marks as broad peaks can have ambiguous peak overlaps. As an alternative, the ENCODE3 pipeline involves a python script overlap_peaks.py for histone marks.

In this tutorial, we will restrict ourselves to finding a simple overlap of peaks in our histone marks sample.

We will use bedtools for this operation. Bedtools is developed by QuinlanLab at Utah University and it assists in genomic arithmetic i.e. intersect, merge, count, complement ,and shuffle genomic intervals from multiple files in different formats such as BAM, BED, GFF/GTF, VCF.

You can download bedtools from https://bedtools.readthedocs.io/en/latest/index.html. After installing and adding it to your path, run the following commands.

```bash

#Find overlap between two replicates
bedtools intersect -a ./macs2/output/H3K4me1_GM12878_rep1_blfilt_peaks.broadPeak -b ./macs2/output/H3K4me1_GM12878_rep2_blfilt_peaks.broadPeak  > ./bedtools/H3K4me1_GM12878_overlaps_blfilt.broadPeak

#Remove non-standard chromosome names
grep -P 'chr[\dXY]+[ \t]' ./bedtools/H3K4me1_GM12878_overlaps_blfilt.broadPeak > ./bedtools/H3K4me1_GM12878_overlaps_filt.broadPeak

```

Bedtools intersect tool evaluates A (file 1) and finds regions that overlap in B (file 2). One thing to note is that if you are working with large files, consider sorting them by chromosome and then by start position, before using bedtools intersect (e.g., sort -k1,1 -k2,2n in.bed > in.sorted.bed for BED files). 
This will lead to a faster and memory efficient computation.


These consensus peak files generated by IDR and BedTools are available in **the tutorial_data folder in the repository.**

--------------------------------------------------------------------

### 8 . Visualization of ChIP-Seq Reads <a name="visualization"></a>

Visualization of data obtained from ChIP-Seq analysis gives us an idea about the data quality, aids in understanding the results, assists in identifying interesting candidate genes, and allows us to integrate the data with other sources of information.

Among the various techniques to visualize enrichment patterns, we will be exploring genomic coverage files and genome browsers using deepTools and IGV.

You can download and install them from https://deeptools.readthedocs.io/en/develop/ and https://software.broadinstitute.org/software/igv/ respectively.

#### 8a. Genome Coverage Files

Visualization of data can be carried out immediately after alignment to reference genome or after peak calling. In the first case, BAM files containing aligned reads can be loaded to genome browsers like the UCSC genome browser.

But in regions with a high density of reads (such as peaks), such a visualization does not scale well.   Instead of visualising every single read, a genomic coverage file can be generated to summarise the number of reads overlapping each position in the genome.

Since sequenced reads represent only ends of the chromatin-immunoprecipitated DNA fragments, reads are extended to 200 bp in 3’ direction as 200 bp is the rough estimate of fragment size selected in ChIP experiments.
This helps smooth the signal and ensures that summits of the peak appear above the binding sites. To compare different samples, read counts at each position are normalised to one million of mapped reads using scaling factors such as  Reads Per Kilobase per Million mapped reads (RPKM), counts per million (CPM), bins per million mapped reads (BPM) and 1x depth (reads per genome coverage, RPGC).

This genomic coverage data can be stored in various formats and here we will be using bigWig format. BigWig format is an indexed binary format that helps display dense, continuous data as a graph or track in the genome browser.

BigWig files are created from wiggle (wig) type files. Wig format contains data elements in equal sized bins and is compressed. Both file types help store coverage data.

To create our bigWig files, we will be using deeptools. deepTools is a suite of Python tools developed for the efficient analysis of high-throughput sequencing data, such as ChIP-seq, RNA-seq etc.

Considering you have added deeptools to your path, run the following command. Remember to edit the ‘-p parameter’ to adapt to the number of cores in your system.

```bash

#For Histone Marks
bamCoverage -b ./alignment/sorted_filtered_bam/H3K4me1_GM12878_rep1_aln_sorted.bam -o ./visualization/bigwig/H3K4me1_GM12878_rep1.bw --binSize 30 --normalizeUsing BPM --smoothLength 300 --extendReads 200 --centerReads -p 4 -bl ./data/blacklist/GRCh38_unified_blacklist.bed

#For Transcription Factor 
bamCoverage -b ./alignment/sorted_filtered_bam/Nfxl1_GM12878_rep1_aln_sorted.bam -o ./visualization/bigwig/Nfxl1_GM12878_rep1.bw --binSize 20 --normalizeUsing BPM --smoothLength 60 --extendReads --centerReads -p 4 -bl ./data/blacklist/GRCh38_unified_blacklist.bed

```

You should run this command for each of your replicate file and not your input file. This will give you six bigWig files: H3K4me1_GM12878_rep1.bw, H3K4me1_GM12878_rep2.bw, Nfxl1_GM12878_rep1.bw, Nfxl1_GM12878_rep2.bw, Nfxl1_K562_rep1.bw, Nfxl1_K562_rep2.bw.

DeepTools offers two commands to create genomic coverage files: bamCoverage and bamCompare. bamCoverage takes a single BAM file and returns a bigWig file. bamCompare normalizes two files to each other (i.e. ChIP sample relative to input) and returns a single bigWig file.

Here, we are using bamCompare only for our ChIP samples and not for input samples. Details of the above-mentioned commands are:
- b/bam: BAM file to process
- o/outFileName: output File name
- binSize/bs: Size of the bins, in bases, for the output of the bigwig file
- normalizeUsing: Possible choices: RPKM, CPM, BPM, RPGC, None
- smoothLength: The smooth length defines a window, larger than the binSize, to average the number of reads.
- extendReads: This parameter allows the extension of reads to fragment size. If set, each read is extended, without exception. Single-end: Requires a user specified value for the final fragment length. Reads that already exceed this fragment length will not be extended. Paired-end: Reads with mates are always extended to match the fragment size defined by the two read mates. Unmated reads, mate reads that map too far apart (>4x fragment length) or even map to different chromosomes are treated like single-end reads. The input of a fragment length value is optional. If no value is specified, it is estimated from the data (mean of the fragment size of all mate reads).
- centerReads: By adding this option, reads are centered with respect to the fragment length. For paired-end data, the read is centered at the fragment length defined by the two ends of the fragment. For single-end data, the given fragment length is used. This option is useful to get a sharper signal around enriched regions.
- p/numberOfProcessors: numberOfProcessors to use.
- bl/ blackListFileName: A BED or GTF file containing regions that should be excluded from all analyses.

#### 8b. Profile Plots and HeatMaps

Profile Plots and Heatmaps give a visual idea about the read density around transcription start sites. Many regulatory elements’ binding sites (that regulate transcription process and switch on/off genes) are located near the transcription sites (TSS) of the target genes.

These plots give a global evaluation of enrichment around the TSS. Before plotting the data, we need an intermediate file that contains data in a format that can be averaged along multiple genomic regions and plotted.

deepTools’ computeMatrix command does just that. It calculates scores per genome regions and stores it in a count matrix. It accepts multiple score files (bigWig format) and multiple regions files (BED format). This tool can also be used to filter and sort regions according to their score.

ComputeMatrix provides two options to calculate values:
- reference-point: calculate values around a reference point such as TSS (eg: +/- 1kb around TSS)
- scaled-regions: all regions in the BED file are stretched or shrunken to the length (in bases) indicated by the user. (eg: values for all genes scaled to 30kb)

We will be using the option reference-point and we will evaluate regions +/- 2.5 kb around TSS.

*** Before executing the following command, download the BED file which contains the coordinates for all genes in the human genome hg38 from UCSC Table Browser (https://genome.ucsc.edu/cgi-bin/hgTables). Store this file as hg38_genes.bed in ./data/chr_bed folder.

Run the following command:

```bash

computeMatrix reference-point --referencePoint TSS -b 2500 -a 2500 -R ./data/chr_bed/hg38_genes.bed -S ./visualization/bigwig/H3K4me1_GM12878_*.bw --skipZeros -o ./visualization/matrix/matrixH3K4me1_TSS_2500.gz -p 4 -bl ./data/blacklist/GRCh38_unified_blacklist.bed

```

This command runs on the pair of replicates for the samples in our dataset, so you would need to run this thrice and in turn you would get three files: matrixH3K4me1_TSS_2500.gz, matrixNfxl1_GM12878_TSS_2500.gz, matrixNfxl1_K562_TSS_2500.gz

Details about the parameters used in this command are:
- referencePoint: Possible choices: TSS, TES, center.
- R/regionsFileName: File name or names, in BED or GTF format, containing the regions to plot.
- S/scoreFileName: bigWig file(s) containing the scores to be plotted.
- skipZeros: Whether regions with only scores of zero should be included or not. Default is to include them.
- o/outFileName: File name to save the gzipped matrix file needed by the “plotHeatmap” and “plotProfile” tools.
- p/numberOfProcessors: Number of processors to use.
- bl/blackListFileName: A BED file containing regions that should be excluded from all analyses.

Next, we will use these count matrices to create profile plot and heatmaps. Please run the following command for each of the matrix.

```bash

plotProfile -m ./visualization/matrix/matrixH3K4me1_TSS_2500.gz -out ./visualization/figures/TSS_H3K4me1_profile.png --perGroup --colors red blue --plotTitle "" --samplesLabel "Rep1" "Rep2" --refPointLabel "TSS" -T "H3K4me1 read density" -z ""
plotHeatmap -m ./visualization/matrix/matrixH3K4me1_TSS_2500.gz  -out ./visualization/figures/TSS_H3K4me1_heatmap.png --whatToShow 'heatmap and colorbar'

```

Here are the profile plots and heatmaps for each of our samples. 

<p> </br> </p>
<p align="center">
<img src="https://github.com/ShrutiBaikerikar/ChipSeq_Tutorial/blob/main/images/visualization_image1.png" width="800" height=400 alt="Profile plot of all samples image"/>
</p>

<p align="center">
     <b>Profile plot of all samples across the genome</b>
</p>

The profile plots demonstrate the signal at the Transcription Start Sites and 2.5 kb upstream and downstream of the TSS where enhancers like H3K4me1 and transcription factors like Nfxl1 are likely to bind.


<p> </br> </p>
<p align="center">
<img src="https://github.com/ShrutiBaikerikar/ChipSeq_Tutorial/blob/main/images/visualization_image2.png" width="800" height=400 alt="Heatmaps of all samples image"/>
</p>

<p align="center">
     <b>Heatmaps of all samples indicating coverage around the TSS across the genome</b>
</p>

The heatmaps show average enrichment for active histone modification H3K4me1 and transcription factor Nfxl1 across all transcription sites in the entire human genome.

#### 8c. Qualitative Assessment using IGV
IGV or Integrative Genomics Viewer is an easy-to-use, interactive tool for the visual exploration of genomic data. It supports flexible integration of all the common types of genomic data and metadata.
It helps compare data obtained from your analysis with previously established genomic annotations.

Download and install IGV from https://software.broadinstitute.org/software/igv/.

To visualize our data in IGV we will need the bigwig files we just generated along with our blacklist filtered broadPeak and narrowPeak file.

Next, follow these steps:
1. Start IGV.
2. Load the Human genome (hg38) into IGV using the dropdown menu at the top left of your screen
3. Load the bigwig files H3K4me1_GM12878_rep1.bw and H3K4me1_Gm12878_rep2.bw
4. Next load the blacklist filtered broadPeak files: H3K4me1_GM12878_rep1_bfilt_peaks.broadPeak and H3K4me1_Gm12878_rep2_bfilt_peaks.broadPeak

Your IGV interface should now look something like the screenshot below. By default, you will be in a zoomed out view.

<p> </br> </p>
<p align="center">
<img src="https://github.com/ShrutiBaikerikar/ChipSeq_tutorial/blob/main/images/visualization_image3.png" width="800" height=400 alt="IGV for H3K4me1_GM12878 image"/>
</p>

<p align="center">
     <b>IGV Visual For H3K4me1_GM12878 .bigwig and .broadPeak files</b>
</p>

The broadPeak and bigWig files demonstrate all the regions where the histone mark binds H3K4me1.

Use the pulldown menu to zoom into chromosome 1. Also remember to Autoscale each track by clicking on the name of each file in the left-hand side panel.

In the box, right next to the one where you select chromosome 1, type the gene name SELL. What do you observe?

<p> </br> </p>
<p align="center">
<img src="https://github.com/ShrutiBaikerikar/ChipSeq_Tutorial/blob/main/images/visualization_image4.png" width="800" height=400 alt="IGV for H3K4me1_GM12878 chromosome 1 image"/>
</p>

<p align="center">
     <b>IGV Visual For H3K4me1_GM12878 sample zoomed on chromosome 1 SELL gene </b>
</p>

SELL gene encodes a cell surface adhesion molecule and the encoded protein contains a C-type lectin-like domain. It is is required for binding and subsequent rolling of leucocytes on endothelial cells, facilitating their migration into secondary lymphoid organs and inflammation sites.

Promoter-enhancer looping is a mechanism by which the chromatin structures form a loop which brings the enhancer (clusters of DNA sequences such as histone modifications which is mostly located away from the genes whose expression it can alter) into close proximity to its target promoter (region of DNA where initiation of transcription takes place).

One example of a cell-type specific loop is anchored at the promoter of the SELL gene that is expressed in GM12878 cell line. It interacts with H3K4me1 and this is observed in the form of peaks near the SELL gene.

5. Finally, let’s visually compare our data to the output from the full dataset from ENCODE . This should be available from the IGV server under the “File” pull-down menu, but on the current system it only shows Annotations as a dataset on the server.
So we will use File -> Load from URL option.

Enter these two URLs one by one to load the bigWig files for the same dataset as available on ENCODE.

https://www.encodeproject.org/files/ENCFF895HSK/@@download/ENCFF895HSK.bigWig

https://www.encodeproject.org/files/ENCFF476LYZ/@@download/ENCFF476LYZ.bigWig

Here is what your screen should look like.

<p> </br> </p>
<p align="center">
<img src="https://github.com/ShrutiBaikerikar/ChipSeq_Tutorial/blob/main/images/visualization_image5.png" width="800" height=400 alt="IGV for H3K4me1_GM12878 with ENCODE data image"/>
</p>

<p align="center">
     <b>IGV Visual For H3K4me1_GM12878 .bigwig and .broadPeak files from our analysis and .bigwig files from ENCODE</b>
</p>

<p> </br> </p>
<p align="center">
<img src="https://github.com/ShrutiBaikerikar/ChipSeq_Tutorial/blob/main/images/visualization_image6.png" width="800" height=400 alt="IGV for H3K4me1_GM12878 chromosome 1 with ENCODE data image"/>
</p>

<p align="center">
     <b>IGV Visual For H3K4me1_GM12878 sample zoomed on chromosome 1 SELL gene for our analysis and ENCODE data</b>
</p>

Follow the same steps as above for transcription factor Nfxl1. Now you will be loading the following files:
* bigWig files: Nfxl1_GM12878_rep1.bw , Nfxl1_Gm12878_rep2.bw , Nfxl1_K562_rep1.bw and Nfxl1_K562_rep2.bw
* narrowPeak files: Nfxl1_GM12878_rep1_bfilt_peaks.narrowPeak , Nfxl1_Gm12878_rep2_bfilt_peaks.narrowPeak , Nfxl1_K562_rep1_bfilt_peaks.narrowPeak and Nfxl1_K562_rep2_bfilt_peaks.narrowPeak

Your IGV interface should now look something like the screenshot below. By default, you will be in a zoomed out view.

<p> </br> </p>
<p align="center">
<img src="https://github.com/ShrutiBaikerikar/ChipSeq_Tutorial/blob/main/images/visualization_image7.png" width="800" height=400 alt="IGV for Nfxl1_GM12878 and Nfxl1_K562 image"/>
</p>

<p align="center">
     <b>IGV Visual For Nfxl1_GM12878 and Nfxl1_K562 .bigwig and .broadPeak files</b>
</p>

Use the pulldown menu to zoom into chromosome 8 and enter gene symbol MYC (MYC Proto-Oncogene). This gene is a proto-oncogene and encodes a nuclear phosphoprotein that plays a role in cell cycle progression, apoptosis and cellular transformation. ChIP signals indicating interaction between NFxl1 and MYC are observed in both cell lines. 

<p> </br> </p>
<p align="center">
<img src="https://github.com/ShrutiBaikerikar/ChipSeq_Tutorial/blob/main/images/visualization_image8.png" width="800" height=400 alt="IGV for Nfxl1_GM12878 and Nfxl1_K562 chromosome 8 image"/>
</p>

<p align="center">
     <b>IGV Visual For Nfxl1_GM12878 and Nfxl1_K562 sample zoomed on chromosome 8 MYC gene </b>
</p>

Further zoom into chromosome 20 and enter gene symbol ZNF217. ChIP signal indicating interaction between NFxl1 and ZNF217 seems to be more enriched in GM12878 than K562 cell line.

<p> </br> </p>
<p align="center">
<img src="https://github.com/ShrutiBaikerikar/ChipSeq_Tutorial/blob/main/images/visualization_image9.png" width="800" height=400 alt="IGV for Nfxl1_GM12878 and Nfxl1_K562 chromosome 20 image"/>
</p>

<p align="center">
     <b>IGV Visual For Nfxl1_GM12878 and Nfxl1_K562 sample zoomed on chromosome 20 ZNF217 gene </b>
</p>

Next move to chromosome 19 and enter gene symbol NFIC (Nuclear Factor I C-type).  The protein encoded by this gene belongs to the CTF/NF-I family. These are dimeric DNA-binding proteins, and function as cellular transcription factors. ChIP signals indicating interaction between NFxl1 and NFIC are observed only in the K562 cell line in this dataset.

<p> </br> </p>
<p align="center">
<img src="https://github.com/ShrutiBaikerikar/ChipSeq_Tutorial/blob/main/images/visualization_image10.png" width="800" height=400 alt="IGV for Nfxl1_GM12878 and Nfxl1_K562 chromosome 19 image"/>
</p>

<p align="center">
     <b>IGV Visual For Nfxl1_GM12878 and Nfxl1_K562 sample zoomed on chromosome 19 NFIC gene </b>
</p>

Finally, we will load the processed files from ENCODE dataset. In the image below tracks for bigWig files for NFXL1_GM12878 are coloured red while NFXL1_K562 are in green and our data is in blue. Use the following URLs:

https://www.encodeproject.org/files/ENCFF112OPN/@@download/ENCFF112OPN.bigWig

https://www.encodeproject.org/files/ENCFF768PGN/@@download/ENCFF768PGN.bigWig

https://www.encodeproject.org/files/ENCFF713FGJ/@@download/ENCFF713FGJ.bigWig

https://www.encodeproject.org/files/ENCFF400AJJ/@@download/ENCFF400AJJ.bigWig


<p> </br> </p>
<p align="center">
<img src="https://github.com/ShrutiBaikerikar/ChipSeq_Tutorial/blob/main/images/visualization_image11.png" width="800" height=400 alt="IGV for Nfxl1_GM12878 and Nfxl1_K562 with ENCODE data image"/>
</p>

<p align="center">
     <b>IGV Visual For Nfxl1_GM12878 and Nfxl1_K562 .bigwig and .broadPeak files from our analysis and .bigwig files from ENCODE</b>
</p>


--------------------------------------------------------------------

### 9. Functional Analysis <a name="functional"></a>

Peak coordinates obtained in our BED files indicate the likely locations of where the proteins of interest (in our case H3K4me1 and Nfxl1) bind to the genome. It is important to study these locations as they can give us an idea as to which genes’ expression is influenced by the histone mark and transcription factor being studied.

#### 9a. Peak Annotation and Genomic Context

Annotation of peaks simply involves assignment of gene symbols and other characteristics that are applicable to the peak coordinates. 

The genomic context of the binding site of the TFs or Histone marks helps understand its function in the cell. Genomic context is generally studied in terms of genomic location and distance from genes.
This make sense, since most cis-regulatory elements are located near the transcription start sites of the target genes.

For this analysis, we will be using R package ChIPSeeker. Apart from annotating genes, it provides visualization functions to assesses peak coverage and profile of peaks.

Open up RStudio and open up the chipseq-project that we created previously. Next, open up a new R script and save it as chipseeker.R This script is available in **the scripts folder in the repository.**

Start by installing the necessary libraries and loading them.

**COMMANDS**

```r
BiocManager::install("ChIPseeker")
BiocManager::install(TxDb.Hsapiens.UCSC.hg38.knownGene)
BiocManager::install(EnsDb.Hsapiens.v86)
BiocManager::install(clusterProfiler)
BiocManager::install(AnnotationDbi)
BiocManager::install(org.Hs.eg.db)

library(ChIPseeker)
library(TxDb.Hsapiens.UCSC.hg38.knownGene)
library(EnsDb.Hsapiens.v86)
library(clusterProfiler)
library(AnnotationDbi)
library(org.Hs.eg.db)

```
Next, we will load our data. We will be using the consensus peaksets that we identified using IDR for Transcription factors Nfxl1_GM12878 and Nfxl1_K562 and using bedtools for histone mark H3K4me1_GM12878.
The files are Nfxl1_GM12878_idr0.05_filt.narrowPeak, Nfxl1_K562_idr0.05_filt.narrowPeak and H3K4me1_GM12878__overlaps_filt.broadPeak.

Copy these files from IDR and bedtools folder to consensus peak folder in chipseq_r folder. Now load these files into the R interface.

```r

#loading consensus peaksets
#Change the file paths as per your requirements

files <- c("chip_seq_analysis/chipseq_r/consensus_peaks/H3K4me1_GM12878_overlaps_filt.broadPeak", "chip_seq_analysis/chipseq_r/consensus_peaks/Nfxl1_GM12878_idr0.05_filt.narrowPeak", "chip_seq_analysis/chipseq_r/consensus_peaks/Nfxl1_K562_idr0.05_filt.narrowPeak")
names(files) <- c("H3K4me1_GM12878", "Nfxl1_GM12878", "Nfxl1_K562")

```

Let’s begin annotating our peaks. We will be using the annotation database EnsDb.Hsapiens.v86. These are generated from ENSEMBL. Since UCSC-style chromosome names are used we have to change the style of the chromosome names from Ensembl to UCSC.

The annotatePeak function performs peak annotation and it uses the nearest gene method to identify nearest Transcription Start Sites (TSS) to the given genomic coordinates. It also allows you to specify a max distance from the TSS to cover the binding site locations accurately.  Here, we specify a distance of 2.5kb upstream and downstream from TSS.

Further, we plot the genomic feature distribution and the distribution of TF-binding loci relative to TSS.

```r

#Getting Annotations for peaks
peakAnno <- lapply(files, annotatePeak, TxDb=edb,tssRegion=c(-2500, 2500), annoDb="org.Hs.eg.db", verbose=FALSE)
peakAnno

#barchart of genomic features for each sample
plotAnnoBar(peakAnno)

#Distribution of HistoneMark/TF-binding loci relative to TSS
plotDistToTSS(peakAnno, title="Distribution of histone mark or transcription factor-binding loci relative to TSS")

```

Here are the plots of genomic feature distribution and Tf-binding loci distribution

<p> </br> </p>
<p align="center">
<img src="https://github.com/ShrutiBaikerikar/ChipSeq_tutorial/blob/main/images/annotations_image1.png" width="800" height=400 alt="Genomic feature distribution image"/>
</p>

<p align="center">
     <b>Genomic Feature Distribution for all samples</b>
</p>


<p> </br> </p>
<p align="center">
<img src="https://github.com/ShrutiBaikerikar/ChipSeq_tutorial/blob/main/images/annotations_image2.png" width="800" height=400 alt="TF-binding loci distribution image"/>
</p>

<p align="center">
     <b>Distribution of TF-binding loci relative to TSS for all samples</b>
</p>


#### 9b. Functional Enrichment

After getting peak annotations, we perform functional enrichment analysis to identify major biological processes or functions influenced by these genes.
For this purpose, we use databases of biological ontologies such as Gene Ontology, KEGG etc.

GO or Gene Ontology provides a description of gene functions for several species. It is organised into three non-overlapping ontologies, which describe a protein’s physiological role (Biological Process), molecular activity (Molecular Function) or position within the cell (Cellular Component).
Based on GO annotations, a list of genes can be tested for enrichment of specific functions. For each GO term, the fraction of genes from the list that are associated with this term is compared to its overall occurrence in the sample dataset to identify terms that are significantly over-represented. Significance is often calculated using a p-value from a hypergeometric test.

Similarly KEGG or Kyoto Encyclopaedia of Genes and Genomes provides pathway-level annotations for genes from various species.

To perform gene enrichment analysis, we will be using R package, clusterProfiler.

First, we begin by retrieving the annotations as a dataframe. We then extract the ENTREZ IDs for each gene in every sample. 

Next we run the enrichment analysis using EnrichGO function where ENTREZ IDs of genes from our samples is compared with the human genome to identify over-represented terms. In this case, we are evaluating only ‘Biological Processes’.
We also use enrichKEGG to perform KEGG pathway enrichment.

Results of these analysis are extracted to dataframes and .csv files. The dotplot function helps us visualize the number of genes associated with the first 30 terms (size) and the gene ratio (# genes related to GO term / total number of sig genes).

```r

#Creating Dataframes of annotation
H3K4me1_GM12878_annot <- data.frame(peakAnno[["H3K4me1_GM12878"]]@anno)
head(H3K4me1_GM12878_annot)

Nfxl1_GM12878_annot <- data.frame(peakAnno[["Nfxl1_GM12878"]]@anno)
head(Nfxl1_GM12878_annot)

Nfxl1_K562_annot <- data.frame(peakAnno[["Nfxl1_K562"]]@anno)
head(Nfxl1_K562_annot)

#Getting EntreZ IDs for each sample
H3K4me1_GM12878_entrez <- H3K4me1_GM12878_annot$ENTREZID
Nfxl1_GM12878_entrez <- Nfxl1_GM12878_annot$ENTREZID
Nfxl1_K562_entrez <- Nfxl1_K562_annot$ENTREZID

#Running Gene-Ontology enrichment analysis
H3K4me1_GM12878_ego <- enrichGO(gene = H3K4me1_GM12878_entrez, 
                keyType = "ENTREZID", 
                OrgDb = org.Hs.eg.db, 
                ont = "BP", 
                pAdjustMethod = "BH", 
                qvalueCutoff = 0.05,
                readable = TRUE)

Nfxl1_GM12878_ego <- enrichGO(gene = Nfxl1_GM12878_entrez, 
                                keyType = "ENTREZID", 
                                OrgDb = org.Hs.eg.db, 
                                ont = "BP", 
                                pAdjustMethod = "BH", 
                                qvalueCutoff = 0.05, 
                                readable = TRUE)

Nfxl1_K562_ego <- enrichGO(gene = Nfxl1_K562_entrez, 
                                keyType = "ENTREZID", 
                                OrgDb = org.Hs.eg.db, 
                                ont = "BP", 
                                pAdjustMethod = "BH", 
                                qvalueCutoff = 0.05, 
                                readable = TRUE)

# Output results from GO analysis to a table
H3K4me1_GM12878_summary <- data.frame(H3K4me1_GM12878_ego)
write.csv(H3K4me1_GM12878_summary, "GOanalysis_H3K4me1_GM12878.csv")

Nfxl1_GM12878_summary <- data.frame(Nfxl1_GM12878_ego)
write.csv(Nfxl1_GM12878_summary, "GOanalysis_Nfxl1_GM12878.csv")

Nfxl1_K562_summary <- data.frame(Nfxl1_K562_ego)
write.csv(Nfxl1_K562_summary, "GOanalysis_Nfxl1_K562.csv")

# Dotplot visualization
dotplot(H3K4me1_GM12878_ego, showCategory=30, title="GO enrichment analysis of H3K4me1_GM12878")
dotplot(Nfxl1_GM12878_ego, showCategory=30, title="GO enrichment analysis of Nfxl1_GM12878")
dotplot(Nfxl1_K562_ego, showCategory=30, title="GO enrichment analysis of Nfxl1_K562")

#KEGG enrichment analysis
H3K4me1_GM12878_ekegg <- enrichKEGG(gene = H3K4me1_GM12878_entrez,
                    organism = 'hsa',
                    pvalueCutoff = 0.05)

Nfxl1_GM12878_ekegg <- enrichKEGG(gene = Nfxl1_GM12878_entrez,
                                    organism = 'hsa',
                                    pvalueCutoff = 0.05)

Nfxl1_K562_ekegg <- enrichKEGG(gene = Nfxl1_K562_entrez,
                                  organism = 'hsa',
                                  pvalueCutoff = 0.05)

# Dotplot visualization
dotplot(H3K4me1_GM12878_ekegg, showCategory=30, title="KEGG enrichment analysis of H3K4me1_GM12878")
dotplot(Nfxl1_GM12878_ekegg, showCategory=30, title="KEGG enrichment analysis of Nfxl1_GM12878")
dotplot(Nfxl1_K562_ekegg, showCategory=30, title="KEGG enrichment analysis of Nfxl1_K562")

```
Here are the various enrichment plots for each sample.

<p> </br> </p>
<p align="center">
<img src="https://github.com/ShrutiBaikerikar/ChipSeq_Tutorial/blob/main/images/annotations_image3.png" width="800" height=400 alt="GO enrichment for H3K4me1_GM12878 image"/>
</p>

<p align="center">
     <b>GO Enrichment Analysis for H3K4me1_GM12878</b>
</p>


<p> </br> </p>
<p align="center">
<img src="https://github.com/ShrutiBaikerikar/ChipSeq_Tutorial/blob/main/images/annotations_image4.png" width="800" height=400 alt="GO enrichment for Nfxl1_GM12878 image"/>
</p>

<p align="center">
     <b>GO Enrichment Analysis for Nfxl1_GM12878</b>
</p>


<p> </br> </p>
<p align="center">
<img src="https://github.com/ShrutiBaikerikar/ChipSeq_Tutorial/blob/main/images/annotations_image5.png" width="800" height=400 alt="GO enrichment for Nfxl1_K562 image"/>
</p>

<p align="center">
     <b>GO Enrichment Analysis for Nfxl1_K562</b>
</p>


<p> </br> </p>
<p align="center">
<img src="https://github.com/ShrutiBaikerikar/ChipSeq_Tutorial/blob/main/images/annotations_image6.png" width="800" height=400 alt="KEGG enrichment for H3K4me1_GM12878 image"/>
</p>

<p align="center">
     <b>KEGG Enrichment Analysis for H3K4me1_GM12878</b>
</p>


<p> </br> </p>
<p align="center">
<img src="https://github.com/ShrutiBaikerikar/ChipSeq_Tutorial/blob/main/images/annotations_image7.png" width="800" height=400 alt="KEGG enrichment for Nfxl1_GM12878 image"/>
</p>

<p align="center">
     <b>KEGG Enrichment Analysis for Nfxl1_GM12878</b>
</p>


<p> </br> </p>
<p align="center">
<img src="https://github.com/ShrutiBaikerikar/ChipSeq_Tutorial/blob/main/images/annotations_image8.png" width="800" height=400 alt="KEGG enrichment for Nfxl1_K562 image"/>
</p>

<p align="center">
     <b>KEGG Enrichment Analysis for Nfxl1_K562</b>
</p>


ClusterProfiler also helps us compare enrichment across samples. We will compare differences in enrichment for our Transcription Factor Nfxl1 in different cell lines GM12878 and K562.

```r

#Comparing enrichment across samples
# Create a list with genes from each sample
genes_compare = lapply(c(peakAnno$Nfxl1_GM12878,peakAnno$Nfxl1_K562), function(i) as.data.frame(i)$ENTREZID)

# Run KEGG analysis
compKEGG <- compareCluster(geneCluster = list("Nfxl1_GM12878" = genes_compare[[1]], "Nfxl1_K562" = genes_compare[[2]]), 
                           fun = "enrichKEGG",
                           organism = "hsa",
                           pvalueCutoff  = 0.05, 
                           pAdjustMethod = "BH")
dotplot(compKEGG, showCategory = 20, title = "KEGG Pathway Enrichment Analysis-Comparison")

```


<p> </br> </p>
<p align="center">
<img src="https://github.com/ShrutiBaikerikar/ChipSeq_Tutorial/blob/main/images/annotations_image9.png" width="800" height=400 alt="KEGG enrichment comparison between Nfxl1_GM12878 and Nfxl1_K562 image"/>
</p>

<p align="center">
     <b>Comparison of KEGG Enrichment Analysis between Nfxl1 in GM12878 cell line and Nfxl1 in K562 cell line</b>
</p>

Similar terms crop up for each cluster and it seems that Nfxl1 impacts pathways related to neutrophil function, alcoholism and lupus in both cell lines but terms related to cell signalling pathways involved in tumor development and immune function seem to more significant in K562 cell line. 

K-562 cells are immortalized lymphoblasts isolated from the bone marrow of a 53-year-old chronic myelogenous leukemia patient and are primarily used to study immune system disorders, cell biology, biochemistry, and erythropoiesis. 

GM12878 is a lymphoblastoid cell line produced from the blood of a female donor with northern and western European ancestry by EBV transformation. This cell line has a relatively normal karyotype. It was selected by the International HapMap Project for whole genome and transcriptome sequencing.

--------------------------------------------------------------------

### 10. Motif Analysis <a name="motif"></a>

Analysing the DNA sequences underlying peak regions or motifs provides insights into the DNA binding preference of the studied transcription factor or histone mark that recurrently bind at neighbouring positions in the genome.

Motif analysis involves two strategies:
• De Novo Motif Discovery: to search without prior assumptions for sequences enriched in peak regions
• Known Motif Search: to scan the peak regions for already defined motifs

In De Novo Motif Discovery, the search is usually performed in a window of 50-200 bp around the TF peak summits or whole regions for histone marks. The searches are either:
- Word-based: all possible k-mers or sequences of length k are exhaustively enumerated to generate consensus motifs that occur with increased frequency in the input sequences.
- Profile-based: sequence alignments are iteratively optimised to obtain the best scoring motifs.

Motifs are present with a high frequency across the genome. Therefore, any enriched motif is always tested against control sequences, either provided by the user or generated by randomisation.

Here we will be using the tool HOMER (Hypergeometric Optimization of Motif EnRichment). It is a collection of command line programs for UNIX-style operating systems written in Perl and C++. 
It contains many useful tools for analyzing ChIP-Seq, GRO-Seq, RNA-Seq, DNase-Seq, Hi-C and numerous other types of functional genomics sequencing data sets.

It takes the genomic coordinates of both target and control regions as input or generates random control regions with the possibility to match the GC content of the target regions. 

After GC matching and eliminating sequence bias, HOMER screens it's library of reliable motifs against the target and background sequences for enrichment, returning motifs enriched with a p-value less than 0.05.
The known motif enrichment is performed first since it is usually faster. Finally, HOMER performs de novo motif discovery search for motifs of length 8, 10, and 12 bp.

We will use findMotifsGenome.pl program which is a wrapper that helps set up the data for analysis using the HOMER motif discovery algorithm.  By default, this will perform de novo motif discovery as well as check the enrichment of known motifs.

Please install HOMER in your system using the following guide. [http://homer.ucsd.edu/homer/introduction/install.html] 
Once you have installed HOMER, also download the necessary packages using configureHomer.pl script. These packages are human-o, human-p and hg38.

Now add HOMER to your path and run the following commands:

```bash

findMotifsGenome.pl ./bedtools/H3K4me1_GM12878_overlaps_filt.broadPeak hg38 ./motif/h3k4me1_gm12878_output/ -size 500 -mask -p 4

findMotifsGenome.pl ./idr/Nfxl1_GM12878_idr0.05_filt.narrowPeak hg38 ./motif/nfxl1_gm12878_output/ -size 200 -mask -p 4

findMotifsGenome.pl ./idr/Nfxl1_GM12878_idr0.05_filt.narrowPeak hg38 ./motif/nfxl1_k562_output/ -size 200 -mask -p 4

```

The parameters used in the above command are:
- Input BED file (narrowPeak or broadPeak files)
- Genome (in our case hg38)
- Output directory
- Size: The size of the region used for motif finding. For TFs, we have used 200 and for histone mark the parameter is set to 500.
- Mask: To mask common repeats are usually in both the target and background sequences which cause biases in the result.
- p: Number of CPU cores to use

Open the homerResults.html file for each sample. Here is what the results would look like for De Novo Motif Discovery.

<p> </br> </p>
<p align="center">
<img src="https://github.com/ShrutiBaikerikar/ChipSeq_Tutorial/blob/main/images/motif_image1.png" width="800" height=400 alt="Motif Analysis for H3K4me1_GM12878 image"/>
</p>

<p align="center">
     <b>De Novo Motif Discovery for H3K4me1_GM12878 sample</b>
</p>


<p> </br> </p>
<p align="center">
<img src="https://github.com/ShrutiBaikerikar/ChipSeq_Tutorial/blob/main/images/motif_image2.png" width="800" height=400 alt="Motif Analysis for Nfxl1_GM12878 image"/>
</p>

<p align="center">
     <b>De Novo Motif Discovery for Nfxl1_GM12878 sample</b>
</p>


<p> </br> </p>
<p align="center">
<img src="https://github.com/ShrutiBaikerikar/ChipSeq_Tutorial/blob/main/images/motif_image2.png" width="800" height=400 alt="Motif Analysis for Nfxl1_K562 image"/>
</p>

<p align="center">
     <b>De Novo Motif Discovery for Nfxl1_K562 sample</b>
</p>


HOMER outputs a ranked list of motifs found in the target sequences. 
For each motif, it indicates the sequence (represented in a logo), a p-value corresponding to the enrichment of this motif in target compared to control sequences as well as the best match for this motif within known motifs.


---------------------------------------------------------------------------------------

## Conclusion <a name="conclusion"></a>

ChIP-sequencing is an experiment based on antibody immunoprecipitation (IP) performed on a population of cells to define protein binding sites in the genome using high-throughput sequencing.
It is used in epigenomic analysis: to study genes prior to their expression.

ChIP-sequencing helps study the activity of regulatory elements such as transcription factors and histone marks which regulate or turn on/off genes. After extracting data about DNA to which the proteins bind from the wet-lab experiement, 
the ChIP-sequencing workflow generally comprises of:
- General quality Control of Reads
- Preprocessing of Reads
- Alignment of Reads to Reference Genome
- Filtering Reads
- Peak calling
- ChIP-specific quality control
- Handling Replicates
- Visualization of ChIP-Seq data
- Downstream Analysis : Differential Binding, Functional analysis, Motif Analysis

In our sample dataset, we studied histone mark H3K4me1 in GM12878 cell line and transcription factor Nfxl1 in GM12878 and K562 cell lines. Given the characterstic bimodal distribution of ChIP-seq reads around the binding site, 
transcription factor Nfxl1 gives rise to sharp peaks while histone mark H3K4me1 gives rise to broad peak as it binds to several regions across the chromosome.

In the functional analysis, we observe that TF Nfxl1 binding loci is mostly within 0-3 kb from the transcription start site (TSS) while for histone mark H3K4me1 we observe that binding sites are at a distance of 5-100 kb from the TSS. 

Pathway analysis suggests that H3K4me1 is involved in biological processes and pathways related to immune system functioning, cancer, cell death, inflammation and cell signalling. Genes regulated by Nfxl1 seem to be involved in chromatin organization and assembly,
T cell activation, bone mineralization, skeletal muscle adaptation, spine development, leucocyte activity and disorders such as colorectal cancer, type II diabetes, Lupus, Alcoholism, Human T cell leukemia virus 1 infection.

The top 3 motifs identified during de novo motif discovery for H3K4me1 reveal its interaction with transcription factors ETS-1, RUNX1 and IRF8 thus confirming its regulatory activity and role as a distal enhancer. The most significant motif identified for Nfxl1 is 
a sequence associated with Zinc Finger Protein 263 (ZNF263) which is associated with chromatin modification and transcriptional repression. 



-------------------------------------------------------------------------------
## Citations <a name="citations_list"></a>

**[1]** Davis CA, Hitz BC, Sloan CA, et al. "**The Encyclopedia of DNA elements (ENCODE): data portal update."** Nucleic Acids Res. 2018;46(D1):D794-D801. doi:10.1093/nar/gkx1081
        [[Research paper](https://www.ncbi.nlm.nih.gov/pmc/articles/PMC5753278/)] Datasets: Bradley Bernstein, Broad, 2011 and Michael Snyder, Stanford, 2016 (https://www.encodeproject.org/)

**[2]** https://github.com/hbctraining/Intro-to-ChIPseq

**[3]** Andrews S. (2010)."**FastQC: a quality control tool for high throughput sequence data.**" [[Source Code](http://www.bioinformatics.babraham.ac.uk/projects/fastqc)]
        
**[4]** TrimGalore, Felix Krueger. [[Source Code](http://www.bioinformatics.babraham.ac.uk/projects/trim_galore)]

**[5]** Danecek P, Bonfield JK, Liddle J, Marshall J, Ohan V, Pollard MO, Whitwham A, Keane T, McCarthy SA, Davies RM, Li H. "**Twelve years of SAMtools and BCFtools**"
        GigaScience (2021) 10(2) giab008 [33590861] PMID: 33590861  [[Research paper](https://pubmed.ncbi.nlm.nih.gov/33590861/)]

**[6]** Artem Tarasov, Albert J. Vilella, Edwin Cuppen, Isaac J. Nijman, Pjotr Prins, "**Sambamba: fast processing of NGS alignment formats**", Bioinformatics, Volume 31, Issue 12, 15 June 2015, Pages 2032–2034, 
        https://doi.org/10.1093/bioinformatics/btv098 [[Research paper](https://academic.oup.com/bioinformatics/article/31/12/2032/214758)]

**[7]** Zhang, Y., Liu, T., Meyer, C.A. et al. "**Model-based Analysis of ChIP-Seq (MACS)**". Genome Biol 9, R137 (2008). https://doi.org/10.1186/gb-2008-9-9-r137
        [[Research paper](https://genomebiology.biomedcentral.com/articles/10.1186/gb-2008-9-9-r137)] [[Source Code](https://github.com/macs3-project/MACS)]

**[8]** Carroll TS, Liang Z, Salama R, Stark R, de Santiago I (in press). “**Impact of artefact removal on ChIP quality metrics in ChIP-seq and ChIP-exo data.**” Frontiers in Genetics. 

**[9]** Stark R, Brown G (2011). "**DiffBind: differential binding analysis of ChIP-Seq peak data.**"  

**[10]** Yu G, Wang L, He Q (2015). “**ChIPseeker: an R/Bioconductor package for ChIP peak annotation, comparison and visualization.**” Bioinformatics, 31(14), 2382-2383. doi: 10.1093/bioinformatics/btv145.

**[11]** Heinz S, Benner C, Spann N, Bertolino E et al. "**Simple Combinations of Lineage-Determining Transcription Factors Prime cis-Regulatory Elements Required for Macrophage and B Cell Identities.**"
        Mol Cell 2010 May 28;38(4):576-589. PMID: 20513432

**[12]** Wu T, Hu E, Xu S, Chen M, Guo P, Dai Z, Feng T, Zhou L, Tang W, Zhan L, Fu x, Liu S, Bo X, Yu G (2021). “**clusterProfiler 4.0: A universal enrichment tool for interpreting omics data.**” 
         The Innovation, 2(3), 100141. doi: 10.1016/j.xinn.2021.100141. 

**[13]** Carlson M (2019). "**org.Hs.eg.db: Genome wide annotation for Human. R package version 3.8.2. **"

**[14]**  Aaron R. Quinlan, Ira M. Hall, "**BEDTools: a flexible suite of utilities for comparing genomic features, Bioinformatics**", Volume 26, Issue 6, 15 March 2010, Pages 841–842, 
         https://doi.org/10.1093/bioinformatics/btq033 [[Research paper](https://academic.oup.com/bioinformatics/article/26/6/841/244688)]

**[15]**  Team BC, Maintainer BP (2019). TxDb.Hsapiens.UCSC.hg38.knownGene: Annotation package for TxDb object(s). R package version 3.4.6. 

**[16]**  Rainer J (2017). EnsDb.Hsapiens.v86: Ensembl based annotation package. R package version 2.99.0. 

**[17]**  Pagès H, Carlson M, Falcon S, Li N (2021). AnnotationDbi: Manipulation of SQLite-based annotations in Bioconductor. R package version 1.56.2, https://bioconductor.org/packages/AnnotationDbi. 

**[18]**  Langmead B, Salzberg SL. Fast gapped-read alignment with Bowtie 2. Nat Methods. 2012;9(4):357-359. Published 2012 Mar 4. doi:10.1038/nmeth.1923 [[Research paper](https://www.ncbi.nlm.nih.gov/pmc/articles/PMC3322381/)]


-----------------------------------------------------------------------------------

## License <a name="license_name"></a>

This project is licensed under the MIT License - see the [LICENSE.md](https://github.com/ShrutiBaikerikar/ChipSeq_tutorial/blob/main/LICENSE) file for details

