# Define treatments ####
initial <- 0.01
template <- c(rep('standard', 24), rep('dfssmt', 8))
dilution <- c(initial, initial/2, initial/4, initial/8,
							initial/16, initial/32, initial/64, initial/128)

# Define frameshift pairs ####
fun.n <- c('n03_fun_fwd-n03_fun_rev.n03_fun_fwd-n03_fun_rev',
						'n04_fun_fwd-n04_fun_rev.n04_fun_fwd-n04_fun_rev',
						'n05_fun_fwd-n05_fun_rev.n05_fun_fwd-n05_fun_rev',
						'n06_fun_fwd-n06_fun_rev.n06_fun_fwd-n06_fun_rev')
gi.n <- c('n06_gi_fwd-n06_gi_rev.n06_gi_fwd-n06_gi_rev',
					'n07_gi_fwd-n07_gi_rev.n07_gi_fwd-n07_gi_rev',
					'n08_gi_fwd-n08_gi_rev.n08_gi_fwd-n08_gi_rev',
					'n09_gi_fwd-n09_gi_rev.n09_gi_fwd-n09_gi_rev')

 # Define index pairs ####
tags <- c('P01-01-A.P01-01-B', 'P01-01-C.P01-01-D',
          'P01-01-E.P01-01-F', 'P01-01-G.P01-01-H') |>
  rep(4) |>
  sort() |>
  strsplit(".", fixed = T) |>
  unlist()

# Randomly sample treatments ####
set.seed(666)
dilution.rand <- sample(dilution) |>
  c(sample(dilution)) |>
  c(sample(dilution)) |>
  c(rep(NA, 8))

fun.rand <- sample(fun.n) |>
	strsplit(".", fixed = T) |>
  unlist() |> 
  c(sample(fun.n) |>
  		strsplit(".", fixed = T) |>
      unlist()) |> c(sample(fun.n) |>
      							 	strsplit(".", fixed = T) |>
                       unlist()) |> c(sample(fun.n) |>
                       							 	strsplit(".", fixed = T) |>
                                        unlist())

gi.rand <- sample(gi.n) |>
	strsplit(".", fixed = T) |>
  unlist() |> 
  c(sample(gi.n) |>
  		strsplit(".", fixed = T) |>
      unlist()) |> c(sample(gi.n) |>
      							 	strsplit(".", fixed = T) |>
                       unlist()) |> c(sample(gi.n) |>
                       							 	strsplit(".", fixed = T) |>
                                        unlist())

# Build the dataframe ####
meta <- data.frame(template,
           tags,
           dilution = dilution.rand,
           fun_n = fun.rand,
           gi_n = gi.rand)

controls <- data.frame(template = rep('control', 4),
                       tags = c('P01-02-A', 'P01-02-B',
                                'P01-02-C', 'P01-02-D'),
                       dilution = rep(NA, 4),
                       fun_n = c('n03_fun_fwd-n03_fun_rev',
                       						'n04_fun_fwd-n04_fun_rev',
                       						'n05_fun_fwd-n05_fun_rev',
                       						'n06_fun_fwd-n06_fun_rev'),
											 gi_n = c('n06_fun_fwd-n06_fun_rev',
											 					'n07_fun_fwd-n07_fun_rev',
											 					'n08_fun_fwd-n08_fun_rev',
											 					'n09_fun_fwd-n09_fun_rev'))

out <- rbind(meta, controls)
out$fun_id <- paste(out$tags, out$fun_n, sep = '-')
out$gi_id <- paste(out$tags, out$gi_n, sep = '-')
out$combo <- paste(out$fun_id, out$gi_id, sep = '.')

write.csv(out, file.path('logs', 'meta.csv'), row.names = F)
read.csv(file.path('logs', 'meta.csv'))
