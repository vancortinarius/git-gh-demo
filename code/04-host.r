# Apply LULU post-clustering curation and assign fungal taxonomy ####

# Accept an argument for the number of threads and check it for validity ####
threads <- commandArgs(T) |> as.integer()

if(is.na(threads) == T){
    stop('Argument was converted to NA')
}
if(length(threads) < 1){
    stop('Please specify the number of threads to launch')
}
if(length(threads) > 1){
    stop('Too many arguments have been provided')
}
if(is.numeric(threads) == F){
    stop('Only numeric arguments are accepted')
}
if(threads < 1){
    stop('At least one thread is needed')
} else {
    cat(threads, 'threads requested', '\n')
}

# Install the lulu github package ####
remotes::install_github('gerverska/lulu@v1.0.0', dependencies = F)

# Load packages ####
library(Biostrings)
library(dada2)
library(lulu)

# Load functions ####
source(file.path('code', '00-functions.r'))

# Create output directories ####
out <- '04-host'
logs <- file.path(out, 'logs')
unlink(out, recursive = T)
dir.create(logs, recursive = T)
system(paste('touch', file.path(out, 'README.md')))
unlink('scratch', recursive = T)
dir.create('scratch')

# Load inputs ####
meta <- read.csv(file.path('data', 'meta.csv'))
fun.tab <- readRDS(file.path('03-denoise', 'fun-seq-tab.rds'))
gi.tab <- readRDS(file.path('03-denoise', 'gi-seq-tab.rds'))

# Define fungal and Gigantea samples ####
# fun.samp <- meta$fun_id |> unique()
# gi.samp <- meta$gi_id |> unique()

# Extract the fungal dataset from the DADA2 output ####
fun.study.tab <- fun.tab#[rownames(fun.tab) %in% fun.samp, ]
fun.study.tab <- fun.study.tab[, colSums(fun.study.tab) > 0]

# Rename OTUs and sequences ####
fun.otus <- paste0('OTU.', 1:ncol(fun.study.tab))
fun.seq <- getSequences(fun.study.tab) |> DNAStringSet()
names(fun.seq) <- fun.otus
colnames(fun.study.tab) <- fun.otus

# Identify ASVs with complete fungal sequences ####
writeXStringSet(fun.seq, file = file.path('scratch', 'itsx.fa'))

itsx.flags <- paste('-i', file.path('scratch', 'itsx.fa'),
                    '-t "tracheophyta"',
                    '--preserve T',
                    '--cpu', threads, 
                    '-o', file.path(logs, 'itsx'),
                    '--only_full T')
system2('ITSx', args = itsx.flags)

# Remove ASVs with incomplete fungal sequences ####
itsx.otus <- readDNAStringSet(file.path(logs, 'itsx.ITS2.fasta')) |> names()
itsx.seq <- fun.seq[names(fun.seq) %in% itsx.otus]
itsx.tab <- fun.study.tab[, names(itsx.seq)]

# Perform LULU post-clustering on the fungal subset ####
fun.lulu <- lulu.clust(tab = itsx.tab, seq = itsx.seq, name = 'fun', min.match = 0.97)

# Predict fungal OTU taxonomy ####
set.seed(666)
taxa <- assignTaxonomy(fun.lulu$seq,
                       file.path('data', 'sh_general_release_dynamic_all_psme-noga_29.11.2022.fasta.gz'), # Potentially update with _s_ version!!!
                       minBoot = 80,
                       tryRC = T,
                       multithread = threads,
                       outputBootstraps = T)

# Set the sequence object as a dataframe ####
fun.lulu$seq <- fun.lulu$seq |> as.data.frame()
colnames(fun.lulu$seq) <- 'sequence'

# Add taxonomic information and metadata to the fungal output ####
tax <- taxa$tax |> as.data.frame()
boot <- taxa$boot |> as.data.frame()
trim.tax <- lapply(tax, gsub, pattern="^[[:alpha:]]__", replacement='') |> data.frame()

rownames(trim.tax) <- rownames(fun.lulu$seq)
rownames(boot) <- rownames(fun.lulu$seq)

fun.lulu$tax <- trim.tax
fun.lulu$boot <- boot
fun.lulu$meta <- meta

# Save the LULU output ####
file.path(out, 'fun-lulu.rds') |> saveRDS(fun.lulu, file = _)

# Remove temporary files ####
unlink('scratch', recursive = T)
