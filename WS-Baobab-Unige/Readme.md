
# Introduction
In this tutorial we will map 1M RNA-seq reads on S. aureus genome with Bowtie on Baobab.



# Prerequisite 

## Download genome

In general, a good place to download a genome sequence together with the corresponding GFF annotation is NCBI Assembly (https://www.ncbi.nlm.nih.gov/assembly/GCF_001281145.1/). For human and mouse, one may however prefer the GENCODE project (https://www.gencodegenes.org). Here we will download *S.aureus* genome from NCBI Assembly database.

    # download the genome of S.aureus, strain SA564 from NCBI Assembly including the GFF annotations
    wget ftp://ftp.ncbi.nlm.nih.gov/genomes/all/GCA/001/281/145/GCA_001281145.1_ASM128114v1/GCA_001281145.1_ASM128114v1_genomic.gff.gz
    wget ftp://ftp.ncbi.nlm.nih.gov/genomes/all/GCA/001/281/145/GCA_001281145.1_ASM128114v1/GCA_001281145.1_ASM128114v1_genomic.fna.gz



## Download RNA-seq reads

Using sratoolkit download from NCBI GEO database 100k RNA-seq reads for *S.aureus* genome and convert into FASTQ format.

    module add UHTS/Analysis/sratoolkit/2.8.0
    fastq-dump --bzip2 --split-3 -v -X 100000 SRR3994405

Alternatively, download the FASTQ files here:
 - [FASTQ for read1](SRR3994405_1.fastq.bz2)
 - [FASTQ for read2](SRR3994405_2.fastq.bz2)



## Map reads to the genome

    module add UHTS/Aligner/bowtie2/2.2.4
    












    
