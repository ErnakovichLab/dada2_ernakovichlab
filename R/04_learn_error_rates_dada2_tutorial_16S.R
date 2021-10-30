#' ### 2. INFER sequence variants
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
#'
#+ include=FALSE
# this is to load in the previous R environment and necessary packages
# if you are running the pipeline in pieces with slurm

# Load DADA2 and required packages
library(dada2); packageVersion("dada2") # the dada2 pipeline
library(ShortRead); packageVersion("ShortRead") # dada2 depends on this
library(dplyr); packageVersion("dplyr") # for manipulating data
library(tidyr); packageVersion("tidyr") # for creating the final graph at the end of the pipeline
library(Hmisc); packageVersion("Hmisc") # for creating the final graph at the end of the pipeline
library(ggplot2); packageVersion("ggplot2") # for creating the final graph at the end of the pipeline
library(plotly); packageVersion("plotly") # enables creation of interactive graphs, especially helpful for quality plots

load(file = "dada2_ernakovich_Renv.RData")
#'
#' In this part of the pipeline dada2 will learn to distinguish error from biological 
#' differences using a subset of our data as a training set. After it understands the 
#' error rates, we will reduce the size of the dataset by combining all identical 
#' sequence reads into "unique sequences". Then, using the dereplicated data and 
#' error rates, dada2 will infer the sequence variants (OTUs) in our data. Finally, 
#' we will merge the coresponding forward and reverse reads to create a list of the 
#' fully denoised sequences and create a sequence table from the result.

#' #### Housekeeping step - set up and verify the file names for the output:
# File parsing
filtFs <- list.files(filtpathF, pattern="fastq.gz", full.names = TRUE)
filtRs <- list.files(filtpathR, pattern="fastq.gz", full.names = TRUE)

# Sample names in order
sample.names <- substring(basename(filtFs), regexpr("_", basename(filtFs)) + 1) # doesn't drop fastq.gz
sample.names <- gsub("_R1_001.fastq.gz", "", sample.names)
sample.namesR <- substring(basename(filtRs), regexpr("_", basename(filtRs)) + 1) # doesn't drop fastq.gz
sample.namesR <- gsub("_R2_001.fastq.gz", "", sample.namesR)

# Double check
if(!identical(sample.names, sample.namesR)) stop("Forward and reverse files do not match.")
names(filtFs) <- sample.names
names(filtRs) <- sample.names

#' #### Learn the error rates
set.seed(100) # set seed to ensure that randomized steps are replicatable

# Learn forward error rates (Notes: randomize default is FALSE)
errF <- learnErrors(filtFs, nbases = 1e8, multithread = TRUE, randomize = TRUE)

# Learn reverse error rates
errR <- learnErrors(filtRs, nbases = 1e8, multithread = TRUE, randomize = TRUE)


#' For NovaSeq Data only!
#' If you have NovaSeq data you'll need to "enforce monotonicity" on the error model.
#' in practice this means taking the errF or errR matrices and assigning a value of Q=40
#' to all entries in the row that are lower than it. 
#' see github issue here for more details: https://github.com/benjjneb/dada2/issues/791
# Define a function to assign the value at Q=40 in the error matrix to all that are 
# lower than it.
make.monotone.decreasing <- function(v) sapply(seq_along(v), function(i) max(v[i:length(v)]))

# Modify the forward and reverse error loess functions
errF.md <- t(apply(getErrors(errF), 1, make.monotone.decreasing))

errR.md <- t(apply(getErrors(errR), 1, make.monotone.decreasing))


# Allow side-by-side plotting of monotonic conversion and original errors
errF_mon <- errF
errF_mon$err_out <- errF.md

errR_mon <- errR
errR_mon$err_out <- errR.md


#' #### Plot Error Rates
#' We want to make sure that the machine learning algorithm is learning the error rates properly. In the plots below, the red line represents what we should expect the learned error rates to look like for each of the 16 possible base transitions (A->A, A->C, A->G, etc.) and the black line and grey dots represent what the observed error rates are. If the black line and the red lines are very far off from each other, it may be a good idea to increase the ```nbases``` parameter. This alows the machine learning algorthim to train on a larger portion of your data and may help imporve the fit.

errF_plot <- plotErrors(errF, nominalQ = TRUE)
errR_plot <- plotErrors(errR, nominalQ = TRUE)

errF_plot
errR_plot

errF_plot_mon <- plotErrors(errF_mon, nominalQ = TRUE)
errR_plot_mon <- plotErrors(errR_mon, nominalQ = TRUE)

#'
# write to disk
saveRDS(errF_plot, paste0(filtpathF, "/errF_plot.rds"))
saveRDS(errR_plot, paste0(filtpathR, "/errR_plot.rds"))

saveRDS(errF_plot_mon, paste0(filtpathF, "/errF_plot_mon.rds"))
saveRDS(errR_plot_mon, paste0(filtpathR, "/errR_plot_mon.rds"))

ggsave(plot = errF_plot, filename = paste0(filtpathF, "/errF_plot.png"), 
       width = 10, height = 10, dpi = "retina")
ggsave(plot = errR_plot, filename = paste0(filtpathF, "/errR_plot.png"), 
       width = 10, height = 10, dpi = "retina")

ggsave(plot = errF_plot_mon, filename = paste0(filtpathF, "/errF_plot_mon.png"),
       width = 10, height = 10, dpi = "retina")
ggsave(plot = errR_plot_mon, filename = paste0(filtpathF, "/errR_plot_mon.png"),
       width = 10, height = 10, dpi = "retina")



#' | <span> |
#' | :--- |
#' | **STOP:** If you are running this on Premise, download the plots generated here (errF_plot.png and errR_plot.png) and verify that the error plots look appropriate. If not, adjust the learnErrors() function and re-run this step with slurm. |
#' | <span> |

#+ include=FALSE
# this is to save the R environment if you are running the pipeline in pieces with slurm
save.image(file = "dada2_ernakovich_Renv.RData")
