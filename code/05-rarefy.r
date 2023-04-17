# Rarefy the data to reduce the influence of read depth on results ####

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
library(ggplot2)

# Load functions ####
source(file.path('code', '00-functions.r'))

# Create output directories ####
out <- '05-rarefy'
logs <- file.path(out, 'logs')
unlink(out, recursive = T)
dir.create(logs, recursive = T)
system(paste('touch', file.path(out, 'README.md')))

# Load and process inputs ####
fun.lulu <- readRDS(file.path('04-compile', 'fun-lulu.rds'))
gi.lulu <- readRDS(file.path('04-compile', 'gi-lulu.rds'))

fun.tab <- fun.lulu$tab
fun.tab$fun_id <- rownames(fun.tab)

gi.tab <- gi.lulu$tab
gi.tab$gi_id <- rownames(gi.tab)

# Combine fungal and Gigantea OTU tables and metadata into a single dataframe ####
full <- merge(fun.lulu$meta, fun.tab, by = 'fun_id', all = T)
full <- merge(full, gi.tab, by = 'gi_id', all = T)
full$combo <- paste(full$tags, full$fun_n, full$gi_n, sep = '-') # Ideally fix some of this upstream!!!
colnames(full) <- colnames(full) |> gsub('-', '', x = _) # Ideally fix this upstream!!!

# Convert NA values post-merge and sum all reads ####
otus <- colnames(full)[grepl('OTU', colnames(full))]

full.tab <- full[, otus]
full.tab[is.na(full.tab) == T] <- 0
full[, otus] <- full.tab
full$reads <- subset(full, select = otus) |> rowSums()

# Examine whether sequencing depth, dilution treatments, and/or frameshift codes are correlated ####
standard <- full |> subset(template == 'standard')
standard$fun_n <- standard$fun_n|> gsub("_fun_fwd|_fun_rev", '', x = _) |> gsub('n', 'N', x = _)
standard$gi_n <- standard$gi_n|> gsub("_gi_fwd|_gi_rev", '', x = _) |> gsub('n', 'N', x = _)

bias <- ggplot(standard, aes(x = log10(dilution), y = log10(reads)
                     )) +
    geom_point(aes(color = fun_n, shape = gi_n), size = 3) +
    scale_x_continuous(n.breaks = 9) +
    scale_y_continuous(n.breaks = 8) +
    xlab('\nLog-transformed NOGA dilution ( prop. NOGA gDNA )') +
    ylab('Log-transformed read depth\n') +
    labs(color = 'Fun frameshift',
         shape = 'Gi frameshift') +
    scale_color_manual(values = color.pal) +
    theme_classic() +
    theme(panel.border = element_blank(),
          panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(),
          text = element_text(size = 14),
          axis.line = element_line(colour = 'black'),
          axis.text.x = element_text(angle = 45, vjust = 0.5),
          axis.title.x = element_text(face = 'bold'),
          axis.title.y = element_text(face = 'bold'))

file.path(logs, 'bias.png') |> ggsave(bias, width = 9, height = 6)
file.path(logs, 'bias.rds') |> saveRDS(bias, file = _)

# Remove samples beneath target sequencing depth prior to rarefaction ####
depth <- 740 # Need to show how we got at this number!!!

filt <- full |> subset(x = _, reads >= depth)
filt.tab <- filt[, otus]
filt.tab$combo <- filt$combo

# Convert the OTU table to long format ####
long <- reshape(filt.tab,
                direction = 'long', varying = otus,
                v.names = 'reads', times = otus,
                timevar = 'OTU', idvar = 'combo')

# Create a row for each read of each OTU in each sample ####
replicas <- lapply(long, rep, long$reads) |> as.data.frame() |>
    subset(x = _, select = -reads)

# Subsample each sample a number of times to a given depth ####
combo <- filt.tab$combo
subsamples <- 1000
set.seed(666)
rare <- lapply(combo, rarefy,
               depth = depth, subsamples = subsamples,
               replicas = replicas) |> do.call(rbind, args = _)

# Create a column for each OTU and remove NA values ####
wide <- reshape(rare,
                direction = 'wide', timevar = 'OTU',
                idvar = c('combo', 'iteration'),
                v.names = 'reads',
                sep = '_')
colnames(wide) <- colnames(wide) |> gsub('reads_', '', x = _)
wide[is.na(wide) == T] <- 0

# Subset the dataframe for OTU columns ####
wide.tab <- wide[, colnames(wide) %in% otus]

# Identify fungal OTUs ####
fun.otus <- colnames(wide.tab)[grepl('fun_OTU', colnames(wide.tab)) == T]

# Calculate statistics for each subsample ####
wide.tab$noga_ra <- wide.tab$fun_OTU.1 / depth
wide.tab$noga_load <- wide.tab$fun_OTU.1 / wide.tab$gi_OTU.1

wide.tab$fun_ra <- rowSums(wide.tab[, fun.otus]) / depth
wide.tab$fun_load <- rowSums(wide.tab[, fun.otus]) / wide.tab$gi_OTU.1

stat.tab <- wide.tab[!(colnames(wide.tab) %in% otus)]
stat.tab$combo <- wide$combo
stat.tab$iteration <- wide$iteration

stats <- colnames(stat.tab)[grepl("ra$|load$", colnames(stat.tab)) == T]

# Lengthen the dataframe to make summarization simpler ####
stat.long <- reshape(stat.tab,
                direction = 'long', varying = stats,
                v.names = 'value', times = stats,
                timevar = 'stat', idvar = c('combo', 'iteration'))

# Summarize each statistic for each sample using the subsamples ####
calc <- by(stat.long,
           INDICES = list(stat.long$combo, stat.long$stat),
           FUN = function(x){data.frame(combo = unique(x$combo),
                                        stat = unique(x$stat),
                                        mean = mean(x$value),
                                        sd = sd(x$value))
               }) |> do.call(rbind, args = _)

# Widen the dataframe to equalize the number of rows and number of samples ####
calc.wide <- reshape(calc,
                     direction = 'wide', timevar = 'stat',
                     idvar = 'combo',
                     sep = '_')

# Combine the rarefied output and the full metadata into a single dataframe ####
meta <- merge(full, calc.wide, by = 'combo', all = T)
tab <- meta |> subset(select = colnames(meta)[grepl('OTU', colnames(meta)) == T])
rownames(tab) <- meta$combo
fun.gi <- list(depth = depth,
               tab = tab,
               meta = meta,
               fun_seq = fun.lulu$seq,
               tax = fun.lulu$tax,
               boot = fun.lulu$boot,
               gi_seq = gi.lulu$seq)

# Write the object out to an .rds file ####
file.path(out, 'fun-gi.rds') |> saveRDS(fun.gi, file = _)
