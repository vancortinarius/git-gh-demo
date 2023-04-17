# Check the validity of the parallelization argument ####
if [[ $# -lt 1 ]]; then
  echo 'Error: Please specify the number of threads to launch'
  exit 1
fi

if [[ $# -gt 1 ]]; then
  echo 'Error: Too many arguments have been provided'
  exit 1
fi

if [[ $1 =~ [^0-9]+ ]]; then
  echo 'Error: Only integer arguments are accepted'
  exit 1
fi

if [[ $1 -lt 1 ]]; then
  echo 'Error: At least one thread is needed'
  exit 1
fi

# Create output directories ####
out='01-demultiplex'
logs="${out}/logs"
rm -r $out
mkdir -p $logs
touch $out/README.md
mkdir -p scratch

# The log file path is relative to the location of the demultiplexed output... ####
pheniqs mux -t $1 -c data/stage2.json -I data/raw -O $out -R logs/pheniqs.txt

# Produce a read count summary for each demultiplexed file ####
fwd=$(find $out -name "*R1.fq.gz" | grep -f data/nano-k.txt)
echo stage2 reads | cat >> $logs/reads.txt

for file in $(echo $fwd); do

    stage2=$(basename $file | awk -F '-' '{print $1"-"$2"-"$3}')
    reads=$(echo $(gzip -cd $file | wc -l) / 4 | bc)

    echo $stage2 $reads | cat >> $logs/reads.txt

done

# Make fastqc reports for each output file ####
fastqc $(find $out -name "*fq.gz") -o scratch -t $1

# Make cumulative multiqc reports for forward and reverse reads ####
multiqc -f -o $logs -n R1.html -i 'Demultiplexed forward reads' scratch/*R1_fastqc.zip
multiqc -f -o $logs -n R2.html -i 'Demultiplexed reverse reads' scratch/*R2_fastqc.zip

# Remove the scratch folder ####
rm -r scratch
