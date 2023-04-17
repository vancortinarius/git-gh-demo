# This script returns the following: ####

# Load packages ####
library(ggplot2)

# Load trimming summary and experiment metadata ####
meta <- read.csv(file.path('data', 'meta.csv')) |>
    pivot_longer(ends_with('id'), values_to = 'id', names_to = 'primers') |>
    select(template, dilution, id, primers, combo)

trim <- read.delim(file.path('logs', '02-trim3-reads.txt'), sep = ' ') |>
		select(id = trim3, reads) |>
		separate(id, c(NA, NA, NA, 'fwd', 'rev'), '-', F) |> 
		mutate(pair = paste0(fwd, '-', rev))

# Identify obvious Fun + Gi chimeras ####
fwd.chimeras <- trim |> filter(grepl('fun', fwd) & grepl('gi', rev)) |> select(reads) |> sum()
rev.chimeras <- trim |> filter(grepl('gi', fwd) & grepl('fun', rev)) |> select(reads) |> sum()
chimeras <- fwd.chimeras + rev.chimeras

# Track shifts in frameshift codes across the entire sequencing run ####
# 5.8S-Fun + ITS4-Fun ####
fun.shifts <- trim |>
		filter(grepl('fun', fwd) & grepl('fun', rev)) |>
		group_by(pair) |> 
		summarize(total = sum(reads)) |>
		arrange(pair)

fun.shifts$pair <- str_remove_all(fun.shifts$pair, "_fun_fwd|_fun_rev") |>
	str_replace_all('n', 'N')

# Calculate the proportion of correct frameshift demultiplexing ####
fun.total <- fun.shifts$total |> sum() + chimeras

fun.pairs <- c('N03-N03', 'N04-N04',
							 'N05-N05', 'N06-N06')
fun.correct <- fun.shifts |>
		filter(pair %in% fun.pairs) |>
		select(total) |>
		sum()

# Calculate rates at which single or double indels occurs in a read ####
fun.del.stop <- which(fun.shifts$pair == 'N03-N03') - 1 

fun.del <- fun.shifts |> slice(1:fun.del.stop)

fun.del.sum <- fun.del |>
		select(total) |>
		sum()

fun.in.start <- which(fun.shifts$pair == 'N06-N06') + 1

fun.in <- fun.shifts |> slice(fun.in.start:n())

fun.in.sum <- fun.in |>
		select(total) |>
		sum()

fun.text <- paste0(fun.correct, ' / ', fun.total, '\nITS2-trimmed reads\ncorrectly demultiplexed', '\n\n',
									 fun.del.sum, ' reads with deletions in N03-N03', '\n',
									 fun.in.sum, ' reads with insertions in N06-N06', '\n\n',
									 'Red and blue dashed lines indicate pairs\nused to calculate deletion and insertion\noccurrence rates, respectively')

# Prepare values for plotting ####
fun.y <- fun.shifts$total |> max() * 0.9

# Overall summary ####
fun.shift.reads <- ggplot(fun.shifts, aes(x = pair, y = total)) +
		geom_bar(stat = 'identity') +
		geom_vline(xintercept = 'N03-N03', color = 'red', linetype = 'dashed') +
		geom_vline(xintercept = 'N06-N06', color = 'blue', linetype = 'dashed') +
		scale_y_continuous(n.breaks = 8) +
		annotate('text', x = 'N02-N06', y = fun.y, label = fun.text) +
		xlab('\nITS2 frameshift pairs') +
		ylab('Total reads\n') +
		theme_cowplot() +
		theme(axis.text.x = element_text(angle = 90, vjust = 0.5),
					axis.title.x = element_text(face = 'bold'),
					axis.title.y = element_text(face = 'bold'))

# GiF + GiR ####
gi.shifts <- trim |>
		filter(grepl('gi', fwd) & grepl('gi', rev)) |>
		group_by(pair) |> 
		summarize(total = sum(reads)) |>
		arrange(pair)

gi.shifts$pair <- str_remove_all(gi.shifts$pair, "_gi_fwd|_gi_rev") |>
	str_replace_all('n', 'N')

# Calculate the proportion of correct frameshift demultiplexing ####
gi.total <- gi.shifts$total |> sum() + chimeras

gi.pairs <- c('N06-N06', 'N07-N07',
							'N08-N08', 'N09-N09')
gi.correct <- gi.shifts |>
	filter(pair %in% gi.pairs) |>
	select(total) |>
	sum()

# Calculate rates at which at least one insertion or deletion occurs in a read ####
gi.del.stop <- which(gi.shifts$pair == 'N06-N06') - 1 

gi.del <- gi.shifts |> slice(1:gi.del.stop)

gi.del.sum <- gi.del |>
	select(total) |>
	sum()

gi.in.start <- which(gi.shifts$pair == 'N09-N09') + 1

gi.in <- gi.shifts |> slice(gi.in.start:n())

gi.in.sum <- gi.in |>
	select(total) |>
	sum()

gi.text <- paste0(gi.correct, ' / ', gi.total, '\nGi-trimmed reads\ncorrectly demultiplexed', '\n\n',
									 gi.del.sum, ' reads with deletions in N06-N06', '\n',
									 gi.in.sum, ' reads with insertions in N09-N09')

# Calculate global chimera rate ####
chimera.text <- paste0('\n', chimeras, ' / ', (fun.total + gi.total - chimeras), ' reads\nas fungal-plant chimeras')
chimera.gi.text <- paste0(gi.text, '\n', chimera.text)

# Prepare values for plotting #### 
gi.y <- gi.shifts$total |> max() * 0.9

# Overall summary ####
gi.shift.reads <- ggplot(gi.shifts, aes(x = pair, y = total)) +
		geom_bar(stat = 'identity') +
		geom_vline(xintercept = 'N06-N06', color = 'red', linetype = 'dashed') +
		geom_vline(xintercept = 'N09-N09', color = 'blue', linetype = 'dashed') +
		annotate('text', x = 'N05-N07', y = gi.y, label = chimera.gi.text) +
		scale_y_continuous(n.breaks = 10) +
		xlab('\nGigantea frameshift pairs') +
		ylab('') +
		theme_cowplot() +
		theme(axis.text.x = element_text(angle = 90, vjust = 0.5),
					axis.title.x = element_text(face = 'bold'),
					axis.title.y = element_text(face = 'bold'))
