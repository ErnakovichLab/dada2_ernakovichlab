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
errF <- learnErrors(filtFs, nbases = 1e10, multithread = TRUE, randomize = TRUE)

# Learn reverse error rates
errR <- learnErrors(filtRs, nbases = 1e10, multithread = TRUE, randomize = TRUE)

saveRDS(errF, paste0(filtpathF, "/errF.rds"))
saveRDS(errR, paste0(filtpathR, "/errR.rds"))
errF <- readRDS(paste0(filtpathF, "/errF.rds"))
errR <- readRDS(paste0(filtpathR, "/errR.rds"))

#' Novaseq trials:
#' 
#' Trial 1 from JacobRPrice
#' alter loess arguments (weights and span and enforce monotonicity)
loessErrfun_mod <- function(trans) {
  qq <- as.numeric(colnames(trans))
  est <- matrix(0, nrow=0, ncol=length(qq))
  for(nti in c("A","C","G","T")) {
    for(ntj in c("A","C","G","T")) {
      if(nti != ntj) {
        errs <- trans[paste0(nti,"2",ntj),]
        tot <- colSums(trans[paste0(nti,"2",c("A","C","G","T")),])
        rlogp <- log10((errs+1)/tot)  # 1 psuedocount for each err, but if tot=0 will give NA
        rlogp[is.infinite(rlogp)] <- NA
        df <- data.frame(q=qq, errs=errs, tot=tot, rlogp=rlogp)
        
        # original
        # ###! mod.lo <- loess(rlogp ~ q, df, weights=errs) ###!
        # mod.lo <- loess(rlogp ~ q, df, weights=tot) ###!
        # #        mod.lo <- loess(rlogp ~ q, df)
        
        # Gulliem Salazar's solution
        # https://github.com/benjjneb/dada2/issues/938
        mod.lo <- loess(rlogp ~ q, df, weights = log10(tot),span = 2)
        
        pred <- predict(mod.lo, qq)
        maxrli <- max(which(!is.na(pred)))
        minrli <- min(which(!is.na(pred)))
        pred[seq_along(pred)>maxrli] <- pred[[maxrli]]
        pred[seq_along(pred)<minrli] <- pred[[minrli]]
        est <- rbind(est, 10^pred)
      } # if(nti != ntj)
    } # for(ntj in c("A","C","G","T"))
  } # for(nti in c("A","C","G","T"))
  
  # HACKY
  MAX_ERROR_RATE <- 0.25
  MIN_ERROR_RATE <- 1e-7
  est[est>MAX_ERROR_RATE] <- MAX_ERROR_RATE
  est[est<MIN_ERROR_RATE] <- MIN_ERROR_RATE
  
  # enforce monotonicity
  # https://github.com/benjjneb/dada2/issues/791
  estorig <- est
  est <- est %>%
    data.frame() %>%
    mutate_all(funs(case_when(. < X40 ~ X40,
                              . >= X40 ~ .))) %>% as.matrix()
  rownames(est) <- rownames(estorig)
  colnames(est) <- colnames(estorig)
  
  # Expand the err matrix with the self-transition probs
  err <- rbind(1-colSums(est[1:3,]), est[1:3,],
               est[4,], 1-colSums(est[4:6,]), est[5:6,],
               est[7:8,], 1-colSums(est[7:9,]), est[9,],
               est[10:12,], 1-colSums(est[10:12,]))
  rownames(err) <- paste0(rep(c("A","C","G","T"), each=4), "2", c("A","C","G","T"))
  colnames(err) <- colnames(trans)
  # Return
  return(err)
}

# check what this looks like
errF_1 <- learnErrors(
  filtFs,
  multithread = TRUE,
  errorEstimationFunction = loessErrfun_mod,
  verbose = TRUE
)
errR_1 <- learnErrors(
  filtRs,
  multithread = TRUE,
  errorEstimationFunction = loessErrfun_mod,
  verbose = TRUE
)

#' Trial 2 enforce monotonicity only:
#'

loessErrfun_mod <- function(trans) {
  qq <- as.numeric(colnames(trans))
  est <- matrix(0, nrow=0, ncol=length(qq))
  for(nti in c("A","C","G","T")) {
    for(ntj in c("A","C","G","T")) {
      if(nti != ntj) {
        errs <- trans[paste0(nti,"2",ntj),]
        tot <- colSums(trans[paste0(nti,"2",c("A","C","G","T")),])
        rlogp <- log10((errs+1)/tot)  # 1 psuedocount for each err, but if tot=0 will give NA
        rlogp[is.infinite(rlogp)] <- NA
        df <- data.frame(q=qq, errs=errs, tot=tot, rlogp=rlogp)
        
        # original
        # ###! mod.lo <- loess(rlogp ~ q, df, weights=errs) ###!
        mod.lo <- loess(rlogp ~ q, df, weights=tot) ###!
        # #        mod.lo <- loess(rlogp ~ q, df)
        
        # Gulliem Salazar's solution
        # https://github.com/benjjneb/dada2/issues/938
        # mod.lo <- loess(rlogp ~ q, df, weights = log10(tot),span = 2)
        
        pred <- predict(mod.lo, qq)
        maxrli <- max(which(!is.na(pred)))
        minrli <- min(which(!is.na(pred)))
        pred[seq_along(pred)>maxrli] <- pred[[maxrli]]
        pred[seq_along(pred)<minrli] <- pred[[minrli]]
        est <- rbind(est, 10^pred)
      } # if(nti != ntj)
    } # for(ntj in c("A","C","G","T"))
  } # for(nti in c("A","C","G","T"))
  
  # HACKY
  MAX_ERROR_RATE <- 0.25
  MIN_ERROR_RATE <- 1e-7
  est[est>MAX_ERROR_RATE] <- MAX_ERROR_RATE
  est[est<MIN_ERROR_RATE] <- MIN_ERROR_RATE
  
  # enforce monotonicity
  # https://github.com/benjjneb/dada2/issues/791
  estorig <- est
  est <- est %>%
    data.frame() %>%
    mutate_all(funs(case_when(. < X40 ~ X40,
                              . >= X40 ~ .))) %>% as.matrix()
  rownames(est) <- rownames(estorig)
  colnames(est) <- colnames(estorig)
  
  # Expand the err matrix with the self-transition probs
  err <- rbind(1-colSums(est[1:3,]), est[1:3,],
               est[4,], 1-colSums(est[4:6,]), est[5:6,],
               est[7:8,], 1-colSums(est[7:9,]), est[9,],
               est[10:12,], 1-colSums(est[10:12,]))
  rownames(err) <- paste0(rep(c("A","C","G","T"), each=4), "2", c("A","C","G","T"))
  colnames(err) <- colnames(trans)
  # Return
  return(err)
}


# check what this looks like
errF_2 <- learnErrors(
  filtFs,
  multithread = TRUE,
  errorEstimationFunction = loessErrfun_mod,
  verbose = TRUE
)

errR_2 <- learnErrors(
  filtRs,
  multithread = TRUE,
  errorEstimationFunction = loessErrfun_mod,
  verbose = TRUE
)

#' Trial 3 alter loess (weights only) and enforce monotonicity
loessErrfun_mod <- function(trans) {
  qq <- as.numeric(colnames(trans))
  est <- matrix(0, nrow=0, ncol=length(qq))
  for(nti in c("A","C","G","T")) {
    for(ntj in c("A","C","G","T")) {
      if(nti != ntj) {
        errs <- trans[paste0(nti,"2",ntj),]
        tot <- colSums(trans[paste0(nti,"2",c("A","C","G","T")),])
        rlogp <- log10((errs+1)/tot)  # 1 psuedocount for each err, but if tot=0 will give NA
        rlogp[is.infinite(rlogp)] <- NA
        df <- data.frame(q=qq, errs=errs, tot=tot, rlogp=rlogp)
        
        # original
        # ###! mod.lo <- loess(rlogp ~ q, df, weights=errs) ###!
        # mod.lo <- loess(rlogp ~ q, df, weights=tot) ###!
        # #        mod.lo <- loess(rlogp ~ q, df)
        
        # Gulliem Salazar's solution
        # https://github.com/benjjneb/dada2/issues/938
        # mod.lo <- loess(rlogp ~ q, df, weights = log10(tot),span = 2)
        
        # only change the weights
        mod.lo <- loess(rlogp ~ q, df, weights = log10(tot))
        
        pred <- predict(mod.lo, qq)
        maxrli <- max(which(!is.na(pred)))
        minrli <- min(which(!is.na(pred)))
        pred[seq_along(pred)>maxrli] <- pred[[maxrli]]
        pred[seq_along(pred)<minrli] <- pred[[minrli]]
        est <- rbind(est, 10^pred)
      } # if(nti != ntj)
    } # for(ntj in c("A","C","G","T"))
  } # for(nti in c("A","C","G","T"))
  
  # HACKY
  MAX_ERROR_RATE <- 0.25
  MIN_ERROR_RATE <- 1e-7
  est[est>MAX_ERROR_RATE] <- MAX_ERROR_RATE
  est[est<MIN_ERROR_RATE] <- MIN_ERROR_RATE
  
  # enforce monotonicity
  # https://github.com/benjjneb/dada2/issues/791
  estorig <- est
  est <- est %>%
    data.frame() %>%
    mutate_all(funs(case_when(. < X40 ~ X40,
                              . >= X40 ~ .))) %>% as.matrix()
  rownames(est) <- rownames(estorig)
  colnames(est) <- colnames(estorig)
  
  # Expand the err matrix with the self-transition probs
  err <- rbind(1-colSums(est[1:3,]), est[1:3,],
               est[4,], 1-colSums(est[4:6,]), est[5:6,],
               est[7:8,], 1-colSums(est[7:9,]), est[9,],
               est[10:12,], 1-colSums(est[10:12,]))
  rownames(err) <- paste0(rep(c("A","C","G","T"), each=4), "2", c("A","C","G","T"))
  colnames(err) <- colnames(trans)
  # Return
  return(err)
}

# check what this looks like
errF_3 <- learnErrors(
  filtFs,
  multithread = TRUE,
  errorEstimationFunction = loessErrfun_mod,
  verbose = TRUE
)


# check what this looks like
errR_3 <- learnErrors(
  filtRs,
  multithread = TRUE,
  errorEstimationFunction = loessErrfun_mod,
  verbose = TRUE
)

#' ## Plotting from trials:
#' 
#' 
# original
errF_plot <- plotErrors(errF, nominalQ = TRUE)
errR_plot <-plotErrors(errR, nominalQ = TRUE)
ggsave(plot = errF_plot, filename = paste0(filtpathF, "/errF_plot.png"), 
       width = 10, height = 10, dpi = "retina")
ggsave(plot = errR_plot, filename = paste0(filtpathF, "/errR_plot.png"), 
       width = 10, height = 10, dpi = "retina")

# Trial 1
errF_plot1 <- plotErrors(errF_1, nominalQ = TRUE)
errR_plot1 <-plotErrors(errR_1, nominalQ = TRUE)
ggsave(plot = errF_plot1, filename = paste0(filtpathF, "/errF_plot1.png"), 
       width = 10, height = 10, dpi = "retina")
ggsave(plot = errR_plot1, filename = paste0(filtpathF, "/errR_plot1.png"), 
       width = 10, height = 10, dpi = "retina")

errF_plot2 <- plotErrors(errF_2, nominalQ = TRUE)
errR_plot2 <-plotErrors(errR_2, nominalQ = TRUE)
ggsave(plot = errF_plot2, filename = paste0(filtpathF, "/errF_plot2.png"), 
       width = 10, height = 10, dpi = "retina")
ggsave(plot = errR_plot2, filename = paste0(filtpathF, "/errR_plot2.png"), 
       width = 10, height = 10, dpi = "retina")

errF_plot3 <- plotErrors(errF_3, nominalQ = TRUE)
errR_plot3 <-plotErrors(errR_3, nominalQ = TRUE)
ggsave(plot = errF_plot3, filename = paste0(filtpathF, "/errF_plot3.png"), 
       width = 10, height = 10, dpi = "retina")
ggsave(plot = errR_plot3, filename = paste0(filtpathF, "/errR_plot3.png"), 
       width = 10, height = 10, dpi = "retina")




#' For NovaSeq Data only!
#' If you have NovaSeq data you'll need to "enforce monotonicity" on the error model.
#' in practice this means taking the errF or errR matrices and assigning a value of Q=40
#' to all entries in the row that are lower than it. 
#' see github issue here for more details: https://github.com/benjjneb/dada2/issues/791
# Define a function to assign the value at Q=40 in the error matrix to all that are 
# lower than it.
# make.monotone.decreasing <- function(v) sapply(seq_along(v), function(i) max(v[i:length(v)]))
# 
# # Modify the forward and reverse error loess functions
# errF.md <- t(apply(getErrors(errF), 1, make.monotone.decreasing))
# 
# errR.md <- t(apply(getErrors(errR), 1, make.monotone.decreasing))
# 
# 
# # Allow side-by-side plotting of monotonic conversion and original errors
# errF_mon <- errF
# errF_mon$err_out <- errF.md
# 
# errR_mon <- errR
# errR_mon$err_out <- errR.md


#' #### Plot Error Rates
#' We want to make sure that the machine learning algorithm is learning the error rates properly. In the plots below, the red line represents what we should expect the learned error rates to look like for each of the 16 possible base transitions (A->A, A->C, A->G, etc.) and the black line and grey dots represent what the observed error rates are. If the black line and the red lines are very far off from each other, it may be a good idea to increase the ```nbases``` parameter. This alows the machine learning algorthim to train on a larger portion of your data and may help imporve the fit.

# errF_plot <- plotErrors(errF, nominalQ = TRUE)
# errR_plot <- plotErrors(errR, nominalQ = TRUE)
# 
# errF_plot
# errR_plot
# 
# errF_plot_mon <- plotErrors(errF_mon, nominalQ = TRUE)
# errR_plot_mon <- plotErrors(errR_mon, nominalQ = TRUE)

#'
# write to disk
# saveRDS(errF_plot, paste0(filtpathF, "/errF_plot.rds"))
# saveRDS(errR_plot, paste0(filtpathR, "/errR_plot.rds"))
# 
# saveRDS(errF_plot_mon, paste0(filtpathF, "/errF_plot_mon.rds"))
# saveRDS(errR_plot_mon, paste0(filtpathR, "/errR_plot_mon.rds"))
# 
# ggsave(plot = errF_plot, filename = paste0(filtpathF, "/errF_plot.png"), 
#        width = 10, height = 10, dpi = "retina")
# ggsave(plot = errR_plot, filename = paste0(filtpathF, "/errR_plot.png"), 
#        width = 10, height = 10, dpi = "retina")
# 
# ggsave(plot = errF_plot_mon, filename = paste0(filtpathF, "/errF_plot_mon.png"),
#        width = 10, height = 10, dpi = "retina")
# ggsave(plot = errR_plot_mon, filename = paste0(filtpathF, "/errR_plot_mon.png"),
#        width = 10, height = 10, dpi = "retina")
# 


#' | <span> |
#' | :--- |
#' | **STOP:** If you are running this on Premise, download the plots generated here (errF_plot.png and errR_plot.png) and verify that the error plots look appropriate. If not, adjust the learnErrors() function and re-run this step with slurm. |
#' | <span> |

#+ include=FALSE
# this is to save the R environment if you are running the pipeline in pieces with slurm
save.image(file = "dada2_ernakovich_Renv.RData")
