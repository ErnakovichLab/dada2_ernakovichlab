#' #### Dereplication, sequence inference, and merging of paired-end reads
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
#' In this part of the pipeline, dada2 will make decisions about assigning sequences to ASVs (called "sequence inference"). There is a major parameter option in the core function dada() that changes how samples are handled during sequence inference. The parameter ```pool = ``` can be set to: ```pool = FALSE``` (default), ```pool = TRUE```, or ```pool = psuedo```. For details on parameter choice, please see below, and further information on this blogpost [http://fiererlab.org/2020/02/17/whats-in-a-number-estimating-microbial-richness-using-dada2/](http://fiererlab.org/2020/02/17/whats-in-a-number-estimating-microbial-richness-using-dada2/), and explanation on the dada2 tutorial [https://benjjneb.github.io/dada2/pool.html](https://benjjneb.github.io/dada2/pool.html).
#' 
#' **Details**   
#' ```pool = FALSE```: Sequence information is not shared between samples. Fast processing time, less sensitivity to rare taxa.   
#' ```pool = psuedo```: Sequence information is shared in a separate "prior" step. Intermediate processing time, intermediate sensitivity to rare taxa.   
#' ```pool = TRUE```: Sequence information from all samples is pooled together. Slow processing time, most sensitivity to rare taxa.   
#' 

#' #### Default: SAMPLES NOT POOLED
#' For simple communities or when you do not need high sensitivity for rare taxa

# make lists to hold the loop output
mergers <- vector("list", length(sample.names))
names(mergers) <- sample.names
ddF <- vector("list", length(sample.names))
names(ddF) <- sample.names
ddR <- vector("list", length(sample.names))
names(ddR) <- sample.names

# For each sample, get a list of merged and denoised sequences
for(sam in sample.names) {
  cat("Processing:", sam, "\n")
  # Dereplicate forward reads
  derepF <- derepFastq(filtFs[[sam]])
  # Infer sequences for forward reads
  dadaF <- dada(derepF, err = errF_4, multithread = TRUE)
  ddF[[sam]] <- dadaF
  # Dereplicate reverse reads
  derepR <- derepFastq(filtRs[[sam]])
  # Infer sequences for reverse reads
  dadaR <- dada(derepR, err = errR_4, multithread = TRUE)
  ddR[[sam]] <- dadaR
  # Merge reads together
  merger <- mergePairs(ddF[[sam]], derepF, ddR[[sam]], derepR)
  mergers[[sam]] <- merger
}

rm(derepF); rm(derepR)


#' #### Alternative: SAMPLES POOLED 
#' For complex communities when you want to preserve rare taxa
#' alternative: swap ```pool = TRUE``` with ```pool = "pseudo"```

#+ eval = FALSE, include=TRUE
# same steps, not in loop

# Dereplicate forward reads
#derepF.p <- derepFastq(filtFs)
#names(derepF.p) <- sample.names
# Infer sequences for forward reads
#dadaF.p <- dada(derepF.p, err = errF, multithread = TRUE, pool = TRUE)
#names(dadaF.p) <- sample.names

# Dereplicate reverse reads
#derepR.p <- derepFastq(filtRs)
#names(derepR.p) <- sample.names
# Infer sequences for reverse reads
#dadaR.p <- dada(derepR.p, err = errR, multithread = TRUE, pool = TRUE)
#names(dadaR.p) <- sample.names

# Merge reads together
#mergers <- mergePairs(dadaF.p, derepF.p, dadaR.p, derepR.p)


#' #### Construct sequence table
#' You will always perform this step whether or not you have pooled or unpooled ASV picking
seqtab <- makeSequenceTable(mergers)

# Save table as an r data object file
dir.create(table.fp)
saveRDS(seqtab, paste0(table.fp, "/seqtab.rds"))

#' | <span> |
#' | :--- |
#' | **STOP - 05_infer_ASVs_dada2_tutorial_16S.R:** If you are running this on Premise, decide if you want the pooled or not-pooled option delete the options you don't want before running this step with slurm. Also make sure to change the error rate model being used if you are not using the default errR and errF. You can change it in the ```dada()``` function option ```err```. Make sure that you change it for both the forward and reverse reads. (You will likely need to change it if you have NovaSeq data.) |
#' | <span> |

#+ include=FALSE
# this is to save the R environment if you are running the pipeline in pieces with slurm
save.image(file = "dada2_ernakovich_Renv.RData")
