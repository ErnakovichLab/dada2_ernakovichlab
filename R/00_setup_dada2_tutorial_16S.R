#' ## Set up (part 2) - You are logged in to premise (or have Rstudio open on your computer) ##
#' First load and test the installed packages to make sure they're working
#'
#' Load DADA2 and required packages
library(dada2); packageVersion("dada2") # the dada2 pipeline
library(ShortRead); packageVersion("ShortRead") # dada2 depends on this
library(dplyr); packageVersion("dplyr") # for manipulating data
library(tidyr); packageVersion("tidyr") # for creating the final graph at the end of the pipeline
library(Hmisc); packageVersion("Hmisc") # for creating the final graph at the end of the pipeline
library(ggplot2); packageVersion("ggplot2") # for creating the final graph at the end of the pipeline
library(plotly); packageVersion("plotly") # enables creation of interactive graphs, especially helpful for quality plots

# Set up pathway to cutadapt (primer trimming tool) and test
cutadapt <- "cutadapt" # CHANGE ME if not on premise; will probably look something like this: "/usr/local/Python27/bin/cutadapt"
system2(cutadapt, args = "--version") # Check by running shell command from R

#' We will now set up the directories for the script. We'll tell the script where our data is, and where we want to put the outputs of the script. We highly recommend NOT putting outputs of this script  directly into your home directory, or into this tutorial directory. A better idea is to create a new project directory to hold the output each project you work on.
#' 

# Set path to shared data folder and contents
data.fp <- "/mnt/home/ernakovich/shared/dada2_tutorial_data"

# List all files in shared folder to check path
list.files(data.fp)

# Set file paths for barcodes file, map file, and fastqs
# Barcodes need to have 'N' on the end of each 12bp sequence for compatability
#map.fp <- file.path(data.fp, "Molecular_Methods_18_515fBC_16S_Mapping_File_SHORT_vFinal_Fierer_10252018.txt")

#' For the tutorial 16S, we will assign taxonomy with Silva db v138, but you might want to use other databases for your data. Below are paths to some of the databases we use often. (If you are on your own computer you can download the database you need from this link [https://benjjneb.github.io/dada2/training.html](https://benjjneb.github.io/dada2/training.html)):
#' 
#'   - 16S bacteria and archaea (SILVA db): ```/mnt/lz01/ernakovich/shared/db_files/dada2/silva_nr99_v138.1_train_set.fa ```
#' 
#'   - ITS fungi (UNITE db): ```/mnt/lz01/ernakovich/shared/db_files/dada2/UNITE_sh_general_release_10.05.2021/sh_general_release_dynamic_10.05.2021.fasta```
#' 
#'   - 18S protists (PR2 db): ```/mnt/lz01/ernakovich/shared/db_files/dada2/pr2_version_4.14.0_SSU_dada2.fasta```
#' 
#' Set file path for the taxonomy database you will use in step 06
db_fp <- "/mnt/lz01/ernakovich/shared/db_files/dada2/silva_nr99_v138.1_train_set.fa" # CHANGE ME, this is silva 138.1, suitable for 16S only

#' Set up file paths in YOUR directory where you want data; 
#' you do not need to create the sub-directories but they are nice to have for organizational purposes. 

project.fp <- "~/dada2_tutorial_test" # CHANGE ME to project directory; don't append with a "/"

dir.create(project.fp)

# Set up names of sub directories to stay organized
preprocess.fp <- file.path(project.fp, "01_preprocess")
filtN.fp <- file.path(preprocess.fp, "filtN")
trimmed.fp <- file.path(preprocess.fp, "trimmed")
filter.fp <- file.path(project.fp, "02_filter") 
table.fp <- file.path(project.fp, "03_tabletax") 

#' | <span> |
#' | :--- |
#' | **STOP - 00_setup_dada2_tutorial_16S.R:** If you are running this on Premise, open up the 00_setup_dada2_tutorial_16S.R script with nano (or your favorite terminal text editor) and adjust the filepaths above appropriately. |
#' | <span> |
#+ include=FALSE
# this is to save the R environment if you are running the pipeline in pieces with slurm, you can safely either run or ignore this line if you are running on your own computer
save.image(file = "dada2_ernakovich_Renv.RData")
