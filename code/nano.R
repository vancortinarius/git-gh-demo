# Define treatments ####
template <- c(rep('standard', 24), rep('dfssmt', 8))
dilution <- c('1', '1/2', '1/4', '1/8', '1/16', '1/32', '1/64', '1/128')
its2.n <- c('N3-N3_N3-N3', 'N4-N4_N4-N4', 'N5-N5_N5-N5', 'N6-N6_N6-N6')
gi.n <- c('N6-N6_N6-N6', 'N7-N7_N7-N7', 'N8-N8_N8-N8', 'N9-N9_N9-N9')

 # Define index pairs ####
tags <- c('P01-01-A_P01-01-B', 'P01-01-C_P01-01-D',
          'P01-01-E_P01-01-F', 'P01-01-G_P01-01-H') |>
  rep(4) |>
  sort() |>
  strsplit('_') |>
  unlist()

# Randomly sample treatments ####
set.seed(666)
dilution.rand <- sample(dilution) |>
  c(sample(dilution)) |>
  c(sample(dilution)) |>
  c(rep(NA, 8))

its2.rand <- sample(its2.n) |>
  strsplit('_') |>
  unlist() |> 
  c(sample(its2.n) |>
      strsplit('_') |>
      unlist()) |> c(sample(its2.n) |>
                       strsplit('_') |>
                       unlist()) |> c(sample(its2.n) |>
                                        strsplit('_') |>
                                        unlist())

gi.rand <- sample(gi.n) |>
  strsplit('_') |>
  unlist() |> 
  c(sample(gi.n) |>
      strsplit('_') |>
      unlist()) |> c(sample(gi.n) |>
                       strsplit('_') |>
                       unlist()) |> c(sample(gi.n) |>
                                        strsplit('_') |>
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
