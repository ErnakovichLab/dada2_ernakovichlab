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
