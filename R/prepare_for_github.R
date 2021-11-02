# Use this file to prepare for github
# Convert R files to rmarkdown
# A note about the "weird" comments in these scripts:
# Throughout these scripts we use Roxygen comments rather than traditional comments. This allows us to have r scripts that can simultaneously be run on slurm and converted (by this script) into Rmarkdown format for display on github as one long tutorial. You can read more about how these work here: [https://bookdown.org/yihui/rmarkdown-cookbook/spin.html](https://bookdown.org/yihui/rmarkdown-cookbook/spin.html). When possible, try not to change them as it may mess with the formatting or the script's ability to be run both on premise and as an r markdown document. Instead use this script to convert the Rscripts into rmarkdown format to make nice reports. 

r_files <- list.files(pattern = '^.*dada2_tutorial_16S.R$', recursive = T)

# Tell the user some info:
writeLines("setting eval to TRUE")
# Make sure that any eval options are TRUE in individual files:
sed <- "sed" 

sed_eval_TRUE <- paste("-i", "'s/knitr::opts_chunk$set(eval = FALSE/knitr::opts_chunk$set(eval = TRUE/g'",
                 r_files)

lapply(sed_eval_TRUE, function(x) {system2(sed, args = x)}) # apply sed command to each of the r_files

writeLines("Creating Rmd files")

rmd_files <- lapply(r_files, function(x) {
  o = knitr::spin(x, knit = FALSE)
  return(o)})
#o = knitr::spin("dada2_tutorial_16S.R", knit = FALSE)

# Create chunks
chunks <- paste0('```{r steps, child = c("', paste(unlist(rmd_files), collapse =  '", "'), '")}\n```\n')

front_matter_knitr_opts <- "

```{r setup, include=FALSE}
# some setup options for outputing markdown files; feel free to ignore these
knitr::opts_chunk$set(eval = TRUE, 
                      include = TRUE, 
                      warning = FALSE, 
                      message = FALSE,
                      collapse = TRUE,
                      dpi = 300,
                      fig.dim = c(9, 9),
                      out.width = '98%',
                      out.height = '98%')

```"
# Wite all steps to a file after knitr opts
write(front_matter_knitr_opts, file = "dada2_tutorial_16S_all.Rmd", append = FALSE)
cat(chunks, sep = '\n', file = "dada2_tutorial_16S_all.Rmd", append = TRUE)


# Knit complete file
writeLines("Rendering file...")
rmarkdown::render("dada2_tutorial_16S_all.Rmd", 
                  output_file = "dada2_tutorial_16S_all.md",
                  output_format = "github_document")

# For convenience change eval back to FALSE when done building
writeLines("setting eval to FALSE")
sed_eval_FALSE <- paste("-i", "'s/knitr::opts_chunk$set(eval = TRUE/knitr::opts_chunk$set(eval = FALSE/g'",
                        r_files)

lapply(sed_eval_FALSE, function(x) {system2(sed, args = x)}) # apply sed command to each of the r_files
