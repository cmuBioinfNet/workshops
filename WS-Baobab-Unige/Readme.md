
# Introduction
In this tutorial we will map 1M RNA-seq reads on S. aureus genome with Bowtie on Baobab.



# Prerequisite 

## Download RNA-seq reads to map

Using sratoolkit download from NCBI GEO database 100k RNA-seq reads for *S.aureus* genome and convert into FASTQ format.

    module add UHTS/Analysis/sratoolkit/2.8.0
    fastq-dump --bzip2 --split-3 -v -X 100000 SRR3994405
    fastq-dump --bzip2 --split-3 -v -X 100000 SRR3994406

Alternatively, download the FASTQ files here:
 - [SRR3994405_1.fastq.bz2](SRR3994405_1.fastq.bz2)
 - [SRR3994405_2.fastq.bz2](SRR3994405_2.fastq.bz2)
 - [SRR3994406_1.fastq.bz2](SRR3994406_1.fastq.bz2)
 - [SRR3994406_2.fastq.bz2](SRR3994406_2.fastq.bz2)


## Download genome

In general, a good place to download a genome sequence together with the corresponding GFF annotation is NCBI Assembly (https://www.ncbi.nlm.nih.gov/assembly/GCF_001281145.1/). For human and mouse, one may however prefer the GENCODE project (https://www.gencodegenes.org). Here we will download *S.aureus* genome from NCBI Assembly database.

    # download the genome of S.aureus, strain MW2 from NCBI Assembly including the GFF annotations
    wget ftp://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/000/011/265/GCF_000011265.1_ASM1126v1/GCF_000011265.1_ASM1126v1_genomic.gff.gz
    wget ftp://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/000/011/265/GCF_000011265.1_ASM1126v1/GCF_000011265.1_ASM1126v1_genomic.fna.gz


## Index the genome with BWA
    module add UHTS/Aligner/bwa/0.7.13
    bwa index GCF_000011265.1_ASM1126v1_genomic.fna.gz
    

# Analysis

## Map the reads to the genome with BWA
    module add UHTS/Aligner/bwa/0.7.13
    module add UHTS/Analysis/samtools/1.3
    bwa mem GCF_000011265.1_ASM1126v1_genomic.fna.gz <(bzcat SRR3994405_1.fastq.bz2) <(bzcat SRR3994405_2.fastq.bz2) | samtools view -Sb - | samtools sort - > SRR3994405.bam
    samtools index SRR3994405.bam
    bwa mem GCF_000011265.1_ASM1126v1_genomic.fna.gz <(bzcat SRR3994406_1.fastq.bz2) <(bzcat SRR3994406_2.fastq.bz2) | samtools view -Sb - | samtools sort - > SRR3994406.bam
    samtools index SRR3994406.bam

## Quantify number of read per gene with R
    module add R/3.3.2;
    R
    
The enter the following R code to quantify the number of read in the genes

    library(GenomicAlignments)
    library(rtracklayer)
    gff <- import.gff("GCF_000011265.1_ASM1126v1_genomic.fna.gz",feature="gene")
    bf <- BamFileList(c("SRR3994405.bam","SRR3994406.bam"))
    n <- summarizeOverlaps(gff,bf)
    
    # list of genes with >=100 mapped reads
    n <- n[rowSums(assay(n)<100)==0]


    

