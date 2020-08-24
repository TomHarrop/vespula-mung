#!/usr/bin/env Rscript

log <- file(snakemake@log[[1]], open = "wt")
sink(log, type = "message")
sink(log, type = "output", append = TRUE)

library(rtracklayer)
library(data.table)

# gff_file <- "data/gff/Vgerm_23_June_2020.gff3"
# error_file <- "output/ms_version/Vgerm/ncbi/ncbi.val"
# my_spec <- "Vgerm"
# fai_file <- "data/fasta/Vgerm_24_August_2020.assembly.annotated.fna.fai"

gff_file <- snakemake@input[["gff"]]
error_file <- snakemake@input[["val"]]
my_spec <- snakemake@wildcards[["spec"]]
fai_file <- snakemake@input[["fai"]]

# read the fai
fai <- fread(fai_file)

# read gff using GenomicRanges
gr <- import.gff3(gff_file)
mc <- as.data.table(mcols(gr, use.names = TRUE))
pc <- mc[, .(Parent = unlist(Parent)), by = ID]
setkey(pc, Parent)

# read the error report
error_lines <- readLines(error_file)
filtered_errors <- error_lines[!grepl("GenCodeMismatch", error_lines)]

# extract the bad records
my_query <- glue::glue(".*({my_spec}s?[[:digit:]]+g[[:digit:]]+).*")
bad_records <- unique(gsub(my_query,
                           "\\1",
                           filtered_errors))


# recurse through MC to get all the childrn
FindAllChildren <- function(x) {
    my_parent <- copy(x)
    children <- c()
    while(pc[my_parent, any(!is.na(ID))]) {
        my_children <- pc[my_parent, ID]
        children <- c(children, my_children)
        my_parent <- my_children
    }
    return(unique(children))
}

bad_children <- unique(unlist(lapply(bad_records, FindAllChildren)))

# remove records not in the fai
rem_from_assembly <- unique(gr[!seqnames(gr) %in% fai[, unique(V1)]]$ID)

# filter the gff on the bad records
bad_ids <- unique(c(bad_records, bad_children, rem_from_assembly))
removed_ids <- bad_ids[bad_ids %in% mc[, ID]]
print(glue::glue("Removing {length(removed_ids)} records"))
print(removed_ids)
keep_ids <- mc[!ID %in% bad_ids, unique(ID)]

x <- gr[gr$ID %in% keep_ids]
export.gff3(x, snakemake@output[["gff"]])

# log
sessionInfo()
