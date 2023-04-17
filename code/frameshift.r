# Functions #####
survey.dist <- function(i){
  
  these <- combos[[i]]
  
  combo <- these |> unlist() |> paste(collapse = ' ')
  
  distances <- diff[row.names(diff) %in% these, colnames(diff) %in% these] |> as.dist()
  
  min_dist <- distances |> min()
  med_dist <- distances |> median()
  mean_dist <- distances |> mean()
  max_dist <- distances |> max()
  
  data.frame(combo = combo,
             min_dist = min_dist,
             med_dist = med_dist,
             max_dist = max_dist)
  
}

# Set parameters of frameshifts ####
shifts <- 6:9
pairs <- expand.grid(fwd = shifts, rev = shifts, stringsAsFactors = F)
row.names(pairs) <- paste('fwd', pairs$fwd, 'rev', pairs$rev, sep = '_')

# Precalculate the distance matrix between frameshifts, preserving integer status ####
diff <- dist(pairs, method = 'manhattan') |> as.matrix()

# Identify all potential 6-partite combinations of frameshift pairs ####
combos <- combn(row.names(pairs), 6, simplify = F)

# For each combination of frameshift pairs, calculate distance summaries ####
compare <- parallel::mclapply(1:length(combos), survey.dist) |> do.call(rbind, args = _)

# Identify the combinations with optimal distance profiles ####
max.compare <- compare |> subset(min_dist == max(min_dist) & med_dist == max(med_dist))
