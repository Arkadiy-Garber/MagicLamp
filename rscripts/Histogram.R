#----------Histogram plot for checking bitscore threshold value
##---------- DO NOT CHANGE ----------##
## Package List
# packages <- c("ggplot2", "dplyr", "tibble", "stringr", "fuzzyjoin", "plotly", "RColorBrewer", "forcats", "argparse")

library.path <- .libPaths()
library("plotly", lib.loc = library.path)
library("ggplot2", lib.loc = library.path)
library("RColorBrewer", lib.loc = library.path)
library("argparse", lib.loc = library.path)
library("dplyr", lib.loc = library.path)
library("forcats", lib.loc = library.path)


## Install packages not yet installed
# installed_packages <- packages %in% rownames(installed.packages())
# if (any(installed_packages == FALSE)) {
#   print(packages[!installed_packages])
#   install.packages(packages[!installed_packages], repos = "http://cran.us.r-project.org")
# }

## Load packages
# invisible(lapply(packages, library, character.only = TRUE))

## function to color code phrases of interest (from "Desc" column)
ff = function(x, patterns, replacements = patterns, fill = NA, ...)
{
  stopifnot(length(patterns) == length(replacements))
  
  ans = rep_len(as.character(fill), length(x))
  empty = seq_along(x)
  
  for(i in seq_along(patterns)) {
    greps = grepl(patterns[[i]], x[empty], ...)
    ans[empty[greps]] = replacements[[i]]
    empty = empty[!greps]
  }
  
  return(ans)
}

## plot theme
plot_theme <- theme(panel.background = element_rect(fill = "white", colour = "black", size = 1, linetype = "solid"),
                    panel.border = element_rect(colour = "black", fill="NA", size=1),
                    panel.grid.major = element_line(size = 0),
                    panel.grid.minor = element_line(size = 0),
                    axis.text = element_text(size=15, colour="black"),
                    axis.title = element_text(face="bold", size=15),
                    legend.position="right",
                    legend.key = element_rect(fill = "white"),
                    legend.title = element_text(face="bold", size=15),
                    legend.text = element_text(size=15))
plot_nomargins_y <- scale_y_continuous(expand = expand_scale(mult = c(0, 0))) # expand_scale = expansion in most recent version of ggplot2...smh
##-----------------------------------##

##---------- CHANGE THE INPUTS ----------##
## set working directory and readin hmmsearch result

args <- commandArgs(trailingOnly = TRUE)

#wd <- args[1]
file <- c(args[2])


##---------------------------------------##

#setwd(wd)
hmm_table <- read.delim(paste0(file), sep="|", stringsAsFactors=FALSE, fileEncoding="latin1")

# hmm_table <- read.delim(file.choose(), sep="|", stringsAsFactors=FALSE, fileEncoding="latin1")

## Make hmm_table numeric and remove any NAs
hmm_table[, c(1:8)] <- sapply(hmm_table[, c(1:8)], as.numeric)
hmm_table <- hmm_table[complete.cases(hmm_table), ]

##---------- CHANGE THE LIST OF INTERESTED PHRASES ----------##
## List of phrases from "Desc" column to color code. The first c() is a list of phrases of interest (partial identification allowed and it is not case-sensitive). The second c() labels these phrases into groups. The label "Other desc" will group all other phrases not specified together (do not delete this). 


# trueKey = c("ssuF", "sulfonate", "molybde")
trueFile = read.csv(args[3], header = FALSE)
trueKey = trueFile$V1
trueVal = rep("TRUE", times=length(trueKey))


# falseKey = c("cell divison", "DNA", "regulator", "spermidine", "TOBE", "integrase")
falseFile = read.csv(args[4], header = FALSE)
falseKey = falseFile$V1
falseVal = rep("FALSE", times=length(falseKey))

hmm_table$color_code <- ff(hmm_table$Desc,
                           c("partial", "hypothetical", "probable", "putative", "unknown", "uncharacterized", trueKey, falseKey),
                           c("partial", "hypothetical", "putative", "putative", "unknown", "unknown", trueVal, falseVal), "Other desc", ignore.case = TRUE)

##-----------------------------------------------------------##

## reorder levels for histogram plot
hmm_table$color_code <- as.factor(hmm_table$color_code)
hmm_table$color_code <- fct_relevel(hmm_table$color_code, "Other desc", after = 0) #put "Other desc" first
hmm_table$color_code <- fct_relevel(hmm_table$color_code, "TRUE", after = Inf) #put "TRUE" last
hmm_table$color_code <- fct_relevel(hmm_table$color_code, "FALSE", after = Inf) #put "FALSE" after "TRUE"

## color scheme
colourCount = length(unique(hmm_table$color_code))
getPalette = colorRampPalette(brewer.pal(8, "Accent"))

## histogram plot (ggplotly creates an interactive plot)

outfile = paste(args[1], "genie-hist.tiff", sep = "/", collapse = NULL)
tiff(outfile, units="in", width=20, height=12, res=150)

hplot <- ggplot(hmm_table, aes(x=Fullseq_score)) +
  geom_histogram(binwidth=1, aes(fill=color_code)) +
  scale_fill_manual(values=getPalette(colourCount)) +
  plot_nomargins_y + plot_theme
hplot
dev.off()

##---------- CAN MANUALLY CHANGE ZOOM-IN LEVEL ----------##
## zoom into histogram plot. Can change ylim and xlim. To view the entire x-range, use xlim = c(0, max(hmm_table$Fullseq_score))

#hplot_zoom <- ggplot(hmm_table, aes(x=Fullseq_score)) +
#  geom_histogram(binwidth=1, aes(fill=color_code)) +
#  coord_cartesian(ylim = c(0,100), xlim = c(0,200)) +
#  scale_fill_manual(values=getPalette(colourCount)) +
#  plot_nomargins_y + plot_theme
#ggplotly(hplot_zoom)

##-------------------------------------------------------##