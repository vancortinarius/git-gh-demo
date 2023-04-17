# Denoise sequences and identify amplicon sequence variants (ASVs) ####

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

# Load packages ####
library(ggplot2) # Should be imported by dada2, but I also thought this about BioStrings...
library(dada2)

# Load functions ####
source(file.path('code', '00-functions.r'))

# Create output directories ####
in.path <- '02-trim'
out <- '03-denoise'
logs <- file.path(out, 'logs')
unlink(out, recursive = T)
dir.create(logs, recursive = T)
dir.create('scratch')
system(paste('touch', file.path(out, 'README.md')))

# Read in forward and reverse reads ####
in.fwd <- in.path |> list.files(pattern = 'R1.fq.gz', full.names = T) |> sort()
in.rev <- in.path |> list.files(pattern = 'R2.fq.gz', full.names = T) |> sort()

# Make paths for trimmed and filtered files ####
filt <- file.path(out, 'filter')
filt.fwd <- gsub(in.path, filt, in.fwd)
filt.rev <- gsub(in.path, filt, in.rev)

# Trim and quality filter ####
trim <- filterAndTrim(in.fwd, filt.fwd,
                      in.rev, filt.rev,
                      maxEE = c(2,2),
                      multithread = threads)

# # Update list of trimmed file paths to exclude samples with no reads passing filters ####
trim.fwd <- filt |> list.files(pattern = 'R1.fq.gz', full.names = T)
trim.rev <- filt |> list.files(pattern = 'R2.fq.gz', full.names = T)

# Make fastqc reports for each sample and read direction after filtering and trimming ####
fastqc.fwd <- paste(paste(trim.fwd, collapse = ' '),
                    '-t', threads,
                    '-o scratch')
fastqc.rev <- paste(paste(trim.rev, collapse = ' '),
                    '-t', threads,
                    '-o scratch')

system2('fastqc', args = fastqc.fwd)
system2('fastqc', args = fastqc.rev)

# Synthesize multiqc reports for each read direction ####
multiqc.fwd <- paste('-f -o', logs,
                     '-n filt-R1.html',
                     '-ip',
                     "-i 'Forward reads after quality filtering and trimming'",
                     'scratch/*R1_fastqc.zip')
multiqc.rev <- gsub('R1.html', 'R2.html', multiqc.fwd) |> 
		gsub('Forward', 'Reverse', x = _) |> 
		gsub('R1_fastqc.zip', 'R2_fastqc.zip', x = _)

system2('multiqc', args = multiqc.fwd)
system2('multiqc', args = multiqc.rev)

# Denoise Fun and Gi reads separately ####
denoise(marker = 'fun')
denoise(marker = 'gi')

# Remove the scratch directory ####
unlink('scratch', recursive = T)
