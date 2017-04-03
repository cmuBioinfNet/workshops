
# Introduction
In this tutorial we will map 1M RNA-seq reads on S. aureus genome with BWA on Baobab cluster.

## Connexion
Connect with SSH using your ISIS login/password. Create a directory for the workshop, and enter in it.

    ssh login@baobab.unige.ch
    mkdir ws
    cd ws

## Documentation
Some documentation on BAOBAB are available [here](http://baobabmaster.unige.ch/enduser/src/enduser/enduser.html)

### Basic SLURM commands

    squeue    # list the running jobs
    srun      # run a command on the cluster
    sbatch    # launch a batch script on the cluster
    skill     # kill a job

### Reservation for the Workshop
5 node of the cluster have been reserved for the workshop. You can use the reserved nodes with the following option:

    srun --partition=shared --reservation=dayera_75 <command>
    
    

## Modules
Similarly to Vital-IT, software versioning is controled with the command `module`.
You can list and search for available modules with `module spider [pattern]`. This command is also useful to list dependencies between modules, then load a module with `module add <module name>`. For example in the following we will use the following modules:

    module add bwa/075a
    module add samtools/1.3
    module add igv
    module add GCC/4.9.3-2.25 OpenMPI/1.10.2 R/3.3.2



# RNA-seq Data Analysis Prerequisite 

## Download RNA-seq reads to map

We want to download 100k RNA-seq reads for *S.aureus* genome using `sratoolkit` in FASTQ format with the following commands:

    fastq-dump --bzip2 --split-3 -v -X 100000 SRR3994405
    fastq-dump --bzip2 --split-3 -v -X 100000 SRR3994406

Unfortunately only sratoolkit version `2.4.5-2` is installed on BAOBAB, and a more recent version is required for the new directory architecture of NCBI to be reconnized. You can install latest `sratoolkit` version in your local folder with:

    curl https://ftp-trace.ncbi.nlm.nih.gov/sra/sdk/current/sratoolkit.current-ubuntu64.tar.gz | tar -zxvf -
    
And then use it:

    ./sratoolkit.2.8.2-1-ubuntu64/bin/fastq-dump --bzip2 --split-3 -v -X 100000 SRR3994405
    ./sratoolkit.2.8.2-1-ubuntu64/bin/fastq-dump --bzip2 --split-3 -v -X 100000 SRR3994406


Alternatively download the FASTQ files from here:
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
    srun bash -c 'bwa index GCF_000011265.1_ASM1126v1_genomic.fna.gz'
    

# Analysis

## Map the reads to the genome with BWA
    srun bash -c 'bwa mem GCF_000011265.1_ASM1126v1_genomic.fna.gz <(bzcat SRR3994405_1.fastq.bz2) <(bzcat SRR3994405_2.fastq.bz2) | samtools view -Sb - | samtools sort - > SRR3994405.bam'
    srun bash -c 'samtools index SRR3994405.bam'
    srun bash -c 'bwa mem GCF_000011265.1_ASM1126v1_genomic.fna.gz <(bzcat SRR3994406_1.fastq.bz2) <(bzcat SRR3994406_2.fastq.bz2) | samtools view -Sb - | samtools sort - > SRR3994406.bam'
    srun bash -c 'samtools index SRR3994406.bam'

## Visualize mapping results 
    # 1) launch X11 server (XQuartz) on your local machine
    # 2) connect with ssh -X to BAOBAB
    ssh -X login@baobab.unige.ch
    srun igv.sh
    
## Quantify number of read per gene with R
    srun R --vanilla
    
### Install few standard Bioconductor packages (in a personal library)
    source("http://bioconductor.org/biocLite.R")
    biocLite(c("GenomicAlignments","rtracklayer"))

### Quantify number of read
Use the following R code to quantify the number of read in the genes

    library(GenomicAlignments)
    library(rtracklayer)
    
    gff <- import.gff("GCF_000011265.1_ASM1126v1_genomic.gff.gz",feature="gene")
    bf <- BamFileList(c("SRR3994405.bam","SRR3994406.bam"))
    n <- summarizeOverlaps(gff,bf)
    
    # uncomment to add a column with the sequence of the start codon
    #values(n)$start_codon <- getSeq(FaFile("GCF_000011265.1_ASM1126v1_genomic.fna"),promoters(rowRanges(n),0,3))
    
    rowRanges(n)$n <- assay(n)
    write.table(as.data.frame(rowRanges(n)),file="count.txt",sep="\t",row.names=FALSE)



    
    


# Useful command

The unix command `disown` is useful to avoid the running jobs you launch in a bash on a server to be killed when you loose the connection.



