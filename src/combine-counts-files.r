#!/usr/bin/Rscript

# Combine counts files 
# Created by: Abdallah Eteleeb <eteleeb@gmail.com>

args = commandArgs(T)
out.dir = args[1]

### combine all counts 
all.cts <- do.call('rbind', lapply(list.files(paste0(out.dir, "/processed_data/counts/"), full.names = TRUE), readRDS))

### write all counts 
saveRDS(all.cts, file=file.path(out.dir, 'processed_data/counts.rds'), compress=FALSE)



