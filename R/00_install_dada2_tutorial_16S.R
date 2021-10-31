#'# dada2 tutorial with NovaSeq dataset for Ernakovich Lab 
#' *This tutorial originally created by Angela Oliverio and Hannah Holland-Moritz. It has been updated for the Ernakovich Lab. Other contributors to this pipeline include: Corinne Walsh, Matt Gebert, and Kunkun Fan*     
#' *Updated October 28, 2021*
#+ include=FALSE
# some setup options for outputing markdown files; feel free to ignore these
knitr::opts_chunk$set(eval = FALSE, 
                      include = TRUE, 
                      warning = FALSE, 
                      message = FALSE,
                      collapse = TRUE,
                      dpi = 300,
                      fig.dim = c(9, 9),
                      out.width = '98%',
                      out.height = '98%')

#' This pipeline runs the dada2 workflow for Big Data (paired-end) with modifications for NovaSeq sequencing base calls
#' 
#' We suggest opening the dada2 tutorial online to understand more about each step. The original pipeline on which this tutorial is based can be found here: [https://benjjneb.github.io/dada2/bigdata_paired.html](https://benjjneb.github.io/dada2/bigdata_paired.html)
#'    
#'
#' | <span> |
#' | :--- |
#' | **NOTE:** There is a slightly different pipeline for ITS and non-"Big data" workflows. The non-"Big data" pipeline, in particular, has very nice detailed explanations for each step and can be found here: [https://benjjneb.github.io/dada2/tutorial.html](https://benjjneb.github.io/dada2/tutorial.html) |
#' | <span> |
#' 
#' ## Preliminary Checklist (part 0) - Before You Begin  ##
#' 
#' 1. Check to make sure you know what your target 'AMPLICON' length. This can vary between primer sets, as well as WITHIN primer sets. For example, ITS (internal transcribed spacer) amplicon can vary from ~100 bps to 300 bps
#'
#'    For examples regarding commonly used primer sets (515f/806r, Fungal ITS2, 1391f/EukBr) see protocols on the Earth Microbiome Project website: [http://press.igsb.anl.gov/earthmicrobiome/protocols-and-standards/](http://press.igsb.anl.gov/earthmicrobiome/protocols-and-standards/)
#' 
#' 2. Check to make sure you know how long your reads should be (i.e., how long should the reads be coming off the sequencer?) This is not the same as fragment length, as many times, especially with longer fragments, the entire fragment
#'    is not being sequenced in one direction. When long _amplicons_ are not sequenced with a _read length_ that allows for substantial overlap between the forward and reverse read, you can potentially insert biases into the data.
#'    If you intend to merge your paired end reads, ensure that your read length is appropriate. For example, with a MiSeq 2 x 150, 300 cycle kit, you will get bidirectional reads of 150 base pairs. 
#'    
#' 3. Make note of which sequencing platform was used, as this can impact both read quality and downstream analysis. In particular, this pipeline is designed to process NovaSeq data which has very different quality scores than HiSeq or MiSeq data.
#' 
#' 4. Decide which database is best suited for your analysis needs. Note that DADA2 requires databases be in a custom format! If a custom database is required, further formatting will be needed to ensure that it can run correctly in dada2.
#'    
#'    See the following link for details regarding database formatting: [https://benjjneb.github.io/dada2/training.html#formatting-custom-databases](https://benjjneb.github.io/dada2/training.html#formatting-custom-databases)  
#' 
#' 5. For additional tutorials and reporting issues, please see link below:    
#'    dada2 tutorial: [https://benjjneb.github.io/dada2/tutorial.html](https://benjjneb.github.io/dada2/tutorial.html)    
#'    dada2 pipeline issues*: [https://github.com/fiererlab/dada2_fiererlab/issues](https://github.com/fiererlab/dada2_fiererlab/issues)     
#'    
#'    *Note by default, only 'OPEN' issues are shown. You can look at all issues by removing "is:open" in the search bar at the top.  
#'    
#' 
#' ## Set up (part 1) - Steps before starting pipeline ##
#' 
#' #### Downloading this tutorial from github
#' Once you have logged in, you can download a copy of the tutorial into your directory on the server. To retrieve the folder with this tutorial from github directly to the server, type the following into your terminal and hit return after each line.
#' 
#' ```bash    
#' wget https://github.com/ernakovichlab/dada2_ernakovichlab/archive/main.zip
#' unzip main.zip
#' ```
#' If there are ever updates to the tutorial on github, you can update the contents of this folder by downloading the new version from the same link as above.
#' 
#' 
#' #### Setup and install software (you will only need to do this the first time, or if you want to update dada2)
#' 
#'    1. Install the conda environment (this will install all the necessary software to run dada2)
#'    2. First start by cleaning up modules, and then loading the anaconda module. 
#'    ```bash
#'    module purge
#'    module load anaconda/colsa
#'    ```
#'    3. Next create a conda local environment that you can use to run the software. This will install everything you need to run dada2.
#'    ```bash
#'    cd dada2_ernakovichlab
#'    conda env create -f dada2_ernakovich.yml
#'    conda activate dada2_ernakovich
#'    ```
#' | <span> |
#' | :--- |
#' | **WARNING:** This installation may take a long time, so only run this code if you have a fairly large chunk of time! |
#' | <span> |
#' 
#' | <span> |
#' | :--- |
#' | **A note about running this on Premise:** To run this on Premise, you will need to submit R-scripts to the job scheduler (slurm). The R scripts in this tutorial can be found in the "R" folder and have been carefully designed so that each step can be run with on slurm with minimal changes. The R scripts are numbered according to their steps. When you are called on to modify a particular step, use a terminal text editor (such as ```nano```) to open up the appropriate R script and edit the code accordingly. For your convenience, there is also a folder called "slurm" which contains ready-made slurm scripts that you can use to submit each R script. The slurm scripts are designed to be submitted from the "slurm" folder. You can submit them by using `cd slurm` to navigate into the slurm folder, and `sbatch xxx_dada2_tutorial_16S.slurm` to submit each script. Throughout this pipeline you will see **STOP** notices. These indicate how you should modify the R script at each stage. |
#' | <span> |
#' 
#' If you are running it on your own computer (runs slower!):
#' 
#' 1. Download this tutorial from github. Go to [the homepage](https://github.com/fiererlab/dada2_fiererlab/dada2_fiererlab), and click the green "Clone or download" button. Then click "Download ZIP", to save it to your computer. Unzip the file to access the R-script.
#' 2. Download the tutorial data from here [http://cme.colorado.edu/projects/bioinformatics-tutorials](http://cme.colorado.edu/projects/bioinformatics-tutorials)
#' 3. Install cutadapt. If you are using conda, you may also use the .yml file to create an environment with cutadapt and all the necessary R packages pre-installed
#'     - cutadapt can be installed from here: [https://cutadapt.readthedocs.io/en/stable/installation.html](https://cutadapt.readthedocs.io/en/stable/installation.html)
#' 4. Download the dada2-formatted reference database of your choice. Link to download here: [https://benjjneb.github.io/dada2/training.html](https://benjjneb.github.io/dada2/training.html)
#'
#' 5. Open the Rmarkdown script in Rstudio. The script is located in the tutorial folder you downloaded in the first step. You can navigate to the proper folder in Rstudio by clicking on the files tab and navigating to the location where you downloaded the github folder. Then click dada2_ernakovichlab and dada2_tutorial_16S_all.Rmd to open the script.
#' 
#' Now, install DADA2 & other necessary packages(if you haven't opted for the conda option). Depending on how you set up Rstudio, you might get a prompt asking if you want to create your own library. Answer 'yes' twice in the console to continue.
#' 
#' | <span> |
#' | :--- |
#' | **WARNING:** This installation may take a long time, so only run this code if these packages are not already installed! |
#' | <span> |
#'
#+ package installation, eval = FALSE, include=TRUE
install.packages("BiocManager")
BiocManager::install("dada2", version = "3.8")

source("https://bioconductor.org/biocLite.R")
biocLite("ShortRead")
install.packages("dplyr")
install.packages("tidyr")
install.packages("Hmisc")
install.packages("ggplot2")
install.packages("plotly")

#' Once the packages are installed, you can check to make sure the auxiliary
#' software is working and set up some of the variables that you will need 
#' along the way.
#'
#' | <span> |
#' | :--- | 
#' | **NOTE:** If you are not working from premise, you will need to change the file paths for cutadapt to where they are stored on your computer/server. |
#' | <span> |
#' 
#' For this tutorial we will be working with some samples that we obtained 16S amplicon data for, from a Illumina Miseq run. The data for these samples can be found on the CME website. [http://cme.colorado.edu/projects/bioinformatics-tutorials](http://cme.colorado.edu/projects/bioinformatics-tutorials)
#' 
