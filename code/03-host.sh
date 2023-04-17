# Find Tracheophyta reads putatively belonging to Pseudotsuga menziesii ####
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

# Set or make directories for input, output and log files ####
in1='03-denoise/filter'
in2='01-demultiplex'
out='03-host'
logs="${out}/logs"
rm -r $out
mkdir -p $logs
touch $out/README.md

# Use ITSxpress to single out Tracheophyta ITS2 (not sure of specificity) ####
for fwd in $(find $in1 -name "*fun_rev-rc-R1.fq.gz"); do

    # Make the name of the corresponding reverse read ####
    rev=$(echo $fwd | sed 's/R1.fq.gz/R2.fq.gz/')
    
    # Make a sample name for each file ####
    sample=$(basename $fwd | awk -F '-' '{print $1"-"$2"-"$3"-"$4"-"$5"-"$6}')

    # Extract ITS ####
    itsxpress --fastq $fwd --fastq2 $rev --region ITS2 --reversed_primers \
    --taxa Tracheophyta --log $logs/${sample}-itsxpress.txt --outfile $out/${sample}.fq.gz --threads $1
    
done
