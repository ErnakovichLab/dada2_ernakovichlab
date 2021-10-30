# Use this file to prepare for github
# Convert R files to rmarkdown
r_files <- list.files(pattern = '^.*dada2_tutorial_16S.R$', recursive = T)

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
rmarkdown::render("dada2_tutorial_16S_all.Rmd", 
                  output_file = "dada2_tutorial_16S_all.md",
                  output_format = "github_document")
