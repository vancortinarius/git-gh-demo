# Colorblind palette ####
color.pal <- c('#000000', '#E69F00', '#56B4E9',
               '#009E73', '#F0E442', '#0072B2',
               '#D55E00', '#CC79A7')

# 03-denoise.r ####
get.n <- function(x){
    sum(getUniques(x))
}

denoise <- function(marker){
    # Set function package requirements ####
    require(dada2)
    
    # Subset files ####
    marker.fwd <- trim.fwd[grepl(marker, trim.fwd) == T]
    marker.rev <- trim.rev[grepl(marker, trim.rev) == T]
    
    # Dereplicate sequences ####
    derep.fwd <- derepFastq(marker.fwd, verbose = F)
    derep.rev <- derepFastq(marker.rev, verbose = F)
    
    # Trim names for each derep object ####
    names(derep.fwd) <- gsub('-rc.R1.fq.gz', '', names(derep.fwd))
    names(derep.rev) <- gsub('-rc.R2.fq.gz', '', names(derep.rev))
    
    # Create downstream references ####
    if(marker == 'fun'){
        insert <- 'Fun'
        
        pairs <- c('N03-N03', 'N04-N04',
                   'N05-N05', 'N06-N06')
        
        stop <- 'N03-N03'
        start <- 'N06-N06'
        
    } else if(marker == 'gi'){
        insert <- 'Gi'
        
        pairs <- c('N06-N06', 'N07-N07',
                   'N08-N08', 'N09-N09')
        
        stop <- 'N06-N06'
        start <- 'N09-N09'
        
    } else {
        stop("Unknown marker specified: please specify 'fun' or 'gi")
    }

    # Learn errors for each read direction ####
    err.fwd <- learnErrors(marker.fwd, multithread = threads, randomize = T)
    err.fwd.plot <- plotErrors(err.fwd, nominalQ = T) + ggtitle(paste(insert, 'forward read error model'))
    file.path(logs, paste0(marker, '-error-fwd.png')) |> ggsave(err.fwd.plot, width = 12, height = 9)
    file.path(logs, paste0(marker, '-error-fwd.rds')) |> saveRDS(err.fwd.plot, file = _)
    err.rev <- learnErrors(marker.rev, multithread = threads, randomize = T)
    err.rev.plot <- plotErrors(err.rev, nominalQ = T) + ggtitle(paste(insert, 'reverse read error model'))
    file.path(logs, paste0(marker, '-error-rev.png')) |> ggsave(err.rev.plot, width = 12, height = 9)
    file.path(logs, paste0(marker, '-error-rev.rds')) |> saveRDS(err.rev.plot, file = _)
    
    # Denoise reads in both directions ####
    dada.fwd <- dada(derep.fwd, err = err.fwd, multithread = threads, pool = 'pseudo')
    dada.rev <- dada(derep.rev, err = err.rev, multithread = threads, pool = 'pseudo')
    
    # Merge forward and reverse reads and make a sequence table ####
    merged <- mergePairs(dada.fwd, derep.fwd,
                         dada.rev, derep.rev,
                         trimOverhang = T)
    seq.tab <- merged |> makeSequenceTable()
    file.path(out, paste0(marker, '-seq-tab.rds')) |> saveRDS(seq.tab, file = _)
    
    # Make a summary log ####
    marker.trim <- trim[grepl(marker, rownames(trim)) == T, ]
    trim.summary <- marker.trim |> data.frame()
    trim.summary$sample <- rownames(trim.summary)
    trim.summary$sample <- gsub('-rc-R1.fq.gz', '', trim.summary$sample)
    
    track <- cbind(sapply(dada.fwd, get.n),
                   sapply(dada.rev, get.n),
                   sapply(merged, get.n)) |>
        data.frame()
    track$sample <- rownames(track)
    
    log <- merge(trim.summary, track, by = 'sample', all.x = T)
    colnames(log) <- c('sample', 'input', 'filtered', 'denoised.fwd', 'denoised.rev', 'merged')
    file.path(logs, paste0(marker, '-reads.txt')) |> write.table(log, file = _, row.names = F, quote = F)
    
    # Calculate sequencing depth ####
    seq.df <- seq.tab |> data.frame()
    
    seq.df$depth <- rowSums(seq.df)
    reads <- seq.df |> subset(select = 'depth')
    
    # Identify frameshift pairs ####
    reads$fwd <- sapply(strsplit(rownames(reads), '-'), '[', 4) |>
        gsub(paste0('_', marker, '_fwd'), '', x = _) |> gsub('n', 'N', x = _)
    reads$rev <- sapply(strsplit(rownames(reads), '-'), '[', 5) |>
        gsub(paste0('_', marker, '_rev'), '', x = _) |> gsub('n', 'N', x = _)
    reads$pair <- paste0(reads$fwd, '-', reads$rev)
    
    # Calculate the read depth of the marker for each frameshift pair ####
    shift <- by(reads,
                INDICES = list(reads$pair),
                FUN = function(x){data.frame(pair = unique(x$pair),
                                             depth = sum(x$depth))
                }) |> do.call(rbind, args = _)
    
    # Calculate the proportion of correct fungal frameshift demultiplexing ####
    total <- shift$depth |> sum()
    correct <- shift |> subset(pair %in% pairs, select = 'depth') |> sum()
    
    del.stop <- which(shift$pair == stop) - 1 
    del.correct <- shift[shift$pair == stop, 'depth'] |> sum()
    del.sum <- shift[1:del.stop, 'depth'] |> sum()
    
    in.start <- which(shift$pair == start) + 1
    in.correct <- shift[shift$pair == start, 'depth'] |> sum()
    in.sum <- shift[in.start:nrow(shift), 'depth'] |> sum()
    
    indel <- data.frame(primers = insert, total = total,
                        correct = correct,
                        del.correct = del.correct, in.correct = in.correct,
                        deletions = del.sum, insertions = in.sum)
    
    file.path(logs, paste0(marker, '-indel.txt')) |> write.table(indel, file = _, row.names = F, quote = F)
    
    # Plot frameshift errors ####
    frameshift <- ggplot(shift, aes(x = pair, y = depth)) +
        geom_bar(stat = 'identity') +
        geom_vline(xintercept = stop, color = 'red', linetype = 'dashed') +
        geom_vline(xintercept = start, color = 'blue', linetype = 'dashed') +
        scale_y_continuous(n.breaks = 8) +
        xlab(paste0('\n', insert, ' frameshift pairs')) +
        ylab('Total reads\n') +
        theme_classic() +
        theme(text = element_text(size = 14),
              axis.text.x = element_text(angle = 90, vjust = 0.5),
              axis.title.x = element_text(face = 'bold'),
              axis.title.y = element_text(face = 'bold'))
    
    file.path(out, paste0(marker, '-frameshift.png')) |> ggsave(frameshift, width = 12, height = 9)
    file.path(out, paste0(marker, '-frameshift.rds')) |> saveRDS(frameshift, file = _)
}

# 04-compile.r ####
lulu.clust <- function(tab, seq, name, min.match = 0.97){
    
    # Collapse identical ASVs ####
    colnames(tab) <- seq |> as.character() |> unname()
    unique.tab <- collapseNoMismatch(tab)
    
    # Remove chimeras ####
    nochim.tab <- removeBimeraDenovo(unique.tab,
                                     method = 'consensus',
                                     multithread = threads,
                                     verbose = T)
    
    # Reassign OTU names ####
    nochim.otus <- paste0(name, '_OTU.', 1:ncol(nochim.tab))
    nochim.seq <- getSequences(nochim.tab) |> DNAStringSet()
    names(nochim.seq) <- nochim.otus
    colnames(nochim.tab) <- nochim.otus
    
    # Prepare data for LULU clustering ####
    t.nochim.tab <- nochim.tab |> t() |> as.data.frame()
    
    # Make a matchlist by aligning all remaining sequences against each other
    writeXStringSet(nochim.seq,
                    file = file.path('scratch', paste0(name, '-vsearch.fa'))
    )
    
    vsearch.flags <- paste('--usearch_global', file.path('scratch', paste0(name, '-vsearch.fa')),
                           '--db', file.path('scratch', paste0(name, '-vsearch.fa')),
                           '--self',
                           '--id', (min.match - 0.01),
                           '--strand plus',
                           '--iddef 1',
                           '--threads', threads,
                           '--userout', file.path('scratch', paste0(name, '-matchlist.txt')),
                           '--log', file.path(logs, paste0(name, '-vsearch.txt')),
                           '--userfields query+target+id',
                           '--maxaccepts 0',
                           '--query_cov 0.9',
                           '--maxhits 10')
    system2('vsearch', args = vsearch.flags)
    
    matchlist <- read.table(file.path('scratch', paste0(name, '-matchlist.txt')),
                            header = F,
                            as.is = T,
                            stringsAsFactors = F)
    
    # Use LULU to make a curation object and update the OTU table and sequences ####
    curation <- lulu(t.nochim.tab, matchlist, minimum_match = (min.match * 100))
    file.path(logs, paste0(name, '-curation.rds')) |> saveRDS(curation, file = _)
    lulu.tab <- curation$curated_table |> t() |> as.data.frame()
    lulu.seq <- nochim.seq[names(nochim.seq) %in% curation$curated_otus]
    
    # Move the LULU log from its original location to the log directory ####
    system(paste('mv lulu.log*', file.path(logs, paste0(name, '-lulu.txt'))
    ))
    
    # Output the new objects ####
    return(list(tab = lulu.tab,
                seq = lulu.seq))
}

# 05-rarefy.r ####
read.count <- function(x, otus, subsample) {
    # Count OTU reads ####
    data.frame(OTU = otus[x],
               reads = length(subsample[subsample %in% otus[x] == T]))
}
resample <- function(x, samp, depth){
    # Randomly sample reads for a sample ####
    subsample <- sample(samp$OTU, size = depth, replace = F)
    
    # Figure out which OTUs are present in the sample and how many there are ####
    otus <- unique(subsample)
    rich <- length(otus)
    
    # Obtain read counts for each OTU ####
    counts <- lapply(1:rich, read.count, otus = otus, subsample = subsample)
    
    # Combine count output into a single dataframe ####
    out <- do.call('rbind', counts)
    out$iteration <- x
    out
}
rarefy <- function(x, depth, subsamples, replicas){
    # Subset reads belonging to a given sample ####
    samp <- subset(replicas, combo == x)
    
    # Sub-sample these reads several times ####
    rare <- lapply(1:subsamples, resample, samp = samp, depth = depth)
    
    # Combine all the output ####
    out <- do.call('rbind', rare)
    out$combo <- x
    out
}
