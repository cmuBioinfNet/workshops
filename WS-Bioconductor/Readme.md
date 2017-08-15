# Introduction
This workshop introduce standard classes of Bioconductor useful for genomic data analysis.

# Installation
In order to have all the same environment, we will use the official Bioconductor docker image 
1. Download and install docker for your platform (<http://docker.com>)
2. Launch the docker application you just install
3. Open settings (or Preferences). In panel *"Advanced"*, set the desired number of CPUs and memory allocated to the VMs
4. Run in a terminal 
   `docker run -p 80:8787 -v $(pwd):/export bioconductor/release_core2:R3.4.1_Bioc3.5`
5. Open your web browser and go the address `http://localhost` (login: rstudio, passwd: rstudio)


# 
  * Get BAM of aligned reads 
  * Get corresponding genome and annotations
  * Determine proportion of the genome covered by a read
  * Make a BigWig of the coverage to import in a genome browser
  

  * `reduce`





