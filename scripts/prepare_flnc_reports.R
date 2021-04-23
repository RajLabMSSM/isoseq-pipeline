
options <- commandArgs(trailingOnly = TRUE)

meta_file <- options[1]

outFile <- options[2]

meta <- readxl::read_excel(meta_file)

library(tidyverse)

all_reports <- 
    map_df(1:nrow(meta), ~{ 
        report_file <- meta$flnc_report_path[.x]
        sample <- meta$sample[.x]
        df <- read_csv(report_file)
        df$primer <- sample
        return(df)
    })

write_csv(all_reports, outFile)

