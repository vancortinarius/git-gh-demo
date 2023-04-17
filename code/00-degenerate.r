# Calculate annealing temperature for each primer degeneracy ####

# Load packages ####
library(ggplot2)

# Load functions ####
source(file.path('code', '00-functions.r'))

# Define an output directory ####
logs <- file.path('data', 'logs')

# Load IDT data ####
fun.degen <- read.csv(file.path('data', 'fun-degen.csv'))

# Plot the distribution of melting temperatures ####
melt.hist <- ggplot(fun.degen, aes(x = melt_c)) +
    geom_histogram(binwidth = 1) +
    scale_x_continuous(n.breaks = 10) +
    scale_y_continuous(n.breaks = 10) +
    xlab('\nMelting temperature (degrees Celsius)') +
    ylab('Counts\n') +
    theme_classic() +
    theme(text = element_text(size = 14),
          axis.title.x = element_text(face = 'bold'),
          axis.title.y = element_text(face = 'bold'),
          legend.title = element_text(face = 'bold'))

file.path(logs, 'melt.png') |> ggsave(melt.hist, width = 9, height = 9)
file.path(logs, 'melt.rds') |> saveRDS(melt.hist, file = _)

degenerate <- function(melt.c = c(57, 59, 61), time.s = 46, diff.c = 3){
    
    df <- fun.degen
    df$melt_c_round <- df$melt_c |> round()
    
    filt <- df |> subset(melt_c_round %in% melt.c)
    
    prof <- by(filt,
               INDICES = filt$melt_c_round,
               FUN = function(x){data.frame(melt_c = unique(x$melt_c_round),
                                            counts = length(x$melt_c_round))
                   }) |> do.call(rbind, args = _)
    
    prof$prop <- prof$counts / sum(prof$counts)
    prof$anneal_c <- prof$melt_c - diff.c
    prof$time_s <- (prof$prop * time.s) |> round()
    
    prof
    
}

# Create output ####
profiles <- degenerate()
file.path(logs, 'profiles.rds') |> saveRDS(profiles, file = _)

