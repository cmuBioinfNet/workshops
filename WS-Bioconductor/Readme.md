# Introduction
This workshop introduces standard classes of Bioconductor useful for genomic data analysis. It tries to show the power of the framework throught the analysis of ChIP-seq data. It proposes to implement a toy peak detection method whose purpose is no to replace state of the art software, but only to show how Bioconductor classes can be used. This workshop was tested under a macOS environment.

# Installation
In order to have all the same environment, we will use the official Bioconductor docker image 
1. Download and install docker for your platform (<http://docker.com>)
2. Launch the docker application you just install
3. Open settings (or Preferences). In panel *"Advanced"*, set the desired number of CPUs and memory allocated to the VMs
4. Run in a terminal 
   `docker run -p 80:8787 -v $(pwd):/export bioconductor/release_core2:R3.4.1_Bioc3.5`.
   This will make the host data of your working directory available in /export of the VM.
5. Open your web browser and go the address `http://localhost` (login: rstudio, passwd: rstudio)


 # Workflow
  ## Data
  * Get BAM of aligned reads 
  * Get corresponding genome and annotations
  * Determine number of aligned read in each BAM
  * Determine proportion of the genome covered by a read
  
  ## Quantify gene expression in RNA-seq data
  
  
  ## Peak detection in ChIP-seq aligned data
  
  ### Peak detection
  * Extend read to 300bp
  * Compute genome coverage of extended reads and normalize
  * Determine regions enriched in read compared to the Control condition
  * Determine regions enriched on both strands
  * Plot cumulative distribution of the peaks size
  * Make a BigWig of the coverage to import in a genome browser
  
  ### Peak annotation
  * retreive sequence of the peaks
  * find the first peak preceeding each gene and determine their distance

  
  
  
  
  





