# Explicitly define every degenerate version of the temperature-limiting 5.8S-Fun primer ####

# Set or make directories for input, output and log files ####
logs="data/logs"
rm -r $logs
mkdir -p $logs

# Define the output directory ####
logs='data/logs'

# Define ambiguous base possibilities ####
Y1='C'
Y2='T'
R1='A'
R2='G'
W1='A'
W2='T'

# Create ambiguity lists ####
Y={$Y1,$Y2}
R={$R1,$R2}
W={$W1,$W2}

# Write a FASTA file ####
eval "echo \>{1..32}" | tr ' ' '\n' > ${logs}/degenerate-names.txt
eval "echo AACTTT"$Y""$R""$R"CAA"$Y"GGATC"$W"CT" | tr ' ' '\n' > ${logs}/degenerate-seq.txt
paste -d '\n' ${logs}/degenerate-names.txt ${logs}/degenerate-seq.txt > ${logs}/degenerate.fa

# Calculate melting temperatures for each degenerate version ####
Rscript code/00-degenerate.r
