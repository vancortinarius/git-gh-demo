# Define treatments ####
template <- c(rep('standard', 24), rep('dfssmt', 8))
dilution <- c(0.01, 0.01/2, 0.01/4, 0.01/8, 0.01/16, 0.01/32, 0.01/64, 0.01/128)
its2.n <- c('n03_fun_fwd-n03_fun_rev.n03_fun_fwd-n03_fun_rev',
						'n04_fun_fwd-n04_fun_rev.n04_fun_fwd-n04_fun_rev',
						'n05_fun_fwd-n05_fun_rev.n05_fun_fwd-n05_fun_rev',
						'n06_fun_fwd-n06_fun_rev.n06_fun_fwd-n06_fun_rev') # change the separating character from '_' to '.'
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

its2.rand <- sample(its2.n) |>
	strsplit(".", fixed = T) |>
  unlist() |> 
  c(sample(its2.n) |>
  		strsplit(".", fixed = T) |>
      unlist()) |> c(sample(its2.n) |>
      							 	strsplit(".", fixed = T) |>
                       unlist()) |> c(sample(its2.n) |>
                       							 	strsplit(".", fixed = T) |>
                                        unlist())

gi.rand <- sample(gi.n) |>
  strsplit('.') |>
  unlist() |> 
  c(sample(gi.n) |>
      strsplit('.') |>
      unlist()) |> c(sample(gi.n) |>
                       strsplit('.') |>
                       unlist()) |> c(sample(gi.n) |>
                                        strsplit('.') |>
                                        unlist())

# Build the dataframe ####
nano <- data.frame(template,
           tags,
           dilution = dilution.rand,
           its2_n = its2.rand,
           gi_n = gi.rand)

controls <- data.frame(template = rep('control', 4),
                       tags = c('P01-02-A', 'P01-02-B',
                                'P01-02-C', 'P01-02-D'),
                       dilution = rep(NA, 4),
                       its2_n = c('N3-N3', 'N4-N4', 'N5-N5', 'N6-N6'),
                       gi_n = c('N6-N6', 'N7-N7', 'N8-N8', 'N9-N9'))

liston <- data.frame(template = rep('liston', 6),
                     tags = c('P01-02-E', 'P01-02-F', 'P01-02-G',
                              'P01-02-H', 'P01-03-A', 'P01-03-B'),
                     dilution = rep(NA, 6),
                     its2_n = rep(NA, 6),
                     gi_n = rep(NA, 6))

out <- rbind(nano, controls, liston)
out$sample <- row.names(out)

stamp <- gsub(' ', '-', gsub(':', '', Sys.time()))

name <- paste0(stamp, '-nano.csv')
write.csv(out, file.path('output', name), row.names = F)
