# Remove 5` gene primers from forward and reverse paired-end sequences. ####
# Trimming is applied with respect to heterogeneity spacer identity.

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
in='01-demultiplex'
out='02-trim'
logs="${out}/logs"
rm -r $out
mkdir -p $logs
touch $out/README.md

# Make directories for scratch data ####
trim1='scratch/trim1'
trim2='scratch/trim2'
trim3='scratch/trim3'
error='scratch/error'
mkdir -p $trim1 $trim2 $trim3 $error

# Define gene primer sequences ####
# Fun ####
fun_fwd='CCTCCGCTTATTGATATGCTTAART' # ITS4-fun
fun_fwd_rc='AYTTAAGCATATCAATAAGCGGAGG' # ITS4-fun, reverse complement
fun_rev='AACTTTYRRCAAYGGATCWCT' #5.8S-fun
fun_rev_rc='AGWGATCCRTTGYYRAAAGTT' # 5.8S-fun, reverse complement
# Gi ####
gi_fwd='CGTGTTCTGGAATATCTACCTC' # GiF
gi_fwd_rc='GAGGTAGATATTCCAGAACACG' # GiF, reverse complement
gi_rev='ATGCCTGAGCACAAACAG' # GiR
gi_rev_rc='CTGTTTGTGCTCAGGCAT' # GiR, reverse complement

# Define frameshift sequences ####
# There might be a simpler way to do this...
n01='N'
n02='NN'
n03='NNN'
n04='NNNN'
n05='NNNNN'
n06='NNNNNN'
n07='NNNNNNN'
n08='NNNNNNNN'
n09='NNNNNNNNN'
n10='NNNNNNNNNN'
n11='NNNNNNNNNNN'

# Create frameshift demultiplexing FASTA files ####
fun_shift_name={n01,n02,n03,n04,n05,n06,n07,n08}
fun_shift_seq={$n01,$n02,$n03,$n04,$n05,$n06,$n07,$n08}
gi_shift_name={n04,n05,n06,n07,n08,n09,n10,n11}
gi_shift_seq={$n04,$n05,$n06,$n07,$n08,$n09,$n10,$n11}

eval "echo \>"$fun_shift_name"_fun_fwd \>"$gi_shift_name"_gi_fwd" | tr ' ' '\n' > scratch/fwd-names.txt
eval "echo "$fun_shift_seq"$fun_fwd "$gi_shift_seq"$gi_fwd" | tr ' ' '\n' > scratch/fwd-seq.txt
paste -d '\n' scratch/fwd-names.txt scratch/fwd-seq.txt > scratch/fwd-adapters.fa

eval "echo \>"$fun_shift_name"_fun_rev \>"$gi_shift_name"_gi_rev" | tr ' ' '\n' > scratch/rev-names.txt
eval "echo "$fun_shift_seq"$fun_rev "$gi_shift_seq"$gi_rev" | tr ' ' '\n' > scratch/rev-seq.txt
paste -d '\n' scratch/rev-names.txt scratch/rev-seq.txt > scratch/rev-adapters.fa

# Demultiplex by forward reads first ####
# This seems like the only way around the "too many open files" error--ideally, we would trim both at the same time
for fwd in $(find $in -name "*R1.fq.gz" | grep -f data/nano-k.txt); do

    # Define the inputs and outputs for the first trim
    rev=$(echo $fwd | sed 's/R1.fq.gz/R2.fq.gz/')
    stage2=$(basename $fwd | awk -F '-' '{print $1"-"$2"-"$3}')

    fwd_trim1="${trim1}/${stage2}-{name}-R1.fq.gz"
    rev_trim1="${trim1}/${stage2}-{name}-R2.fq.gz"

    # json files are empty, so not sure this is working properly--use fastqc until this is fixed
    #json1="${logs}trim1/${stage2}.cutadapt.json"

    cutadapt -e 0.1 --no-indels -j $1 --discard-untrimmed \
    -g ^file:scratch/fwd-adapters.fa \
    -o $fwd_trim1 -p $rev_trim1 $fwd $rev > $logs/trim1-cutadapt.txt

done

# Initiate the read count summary for the first trimming step ####
echo trim1 reads | cat >> $logs/trim1-reads.txt

# Then demultiplex by reverse reads for each forward-demultiplexed file ####
for fwd_trim1_in in $(find $trim1 -name "*fwd-R1.fq.gz"); do

    rev_trim1_in=$(echo $fwd_trim1_in | sed 's/R1.fq.gz/R2.fq.gz/')

    # Count the number of reads in each file after the first trimming step
    if [[ -f $fwd_trim1_in && -f $rev_trim1_in ]]; then
        trim1sum=$(echo $(gzip -cd $fwd_trim1_in | wc -l) / 4 | bc)
    else
        trim1sum=0
    fi

    # Only attempt to demultiplex by reverse reads if a forward-demultiplexed file exists
    if [ $trim1sum -gt 0 ]; then

        echo $fwd_trim1_in | cat >> $trim1/files.txt
        echo $rev_trim1_in | cat >> $trim1/files.txt

        mid=$(basename $fwd_trim1_in | awk -F '-' '{print $1"-"$2"-"$3"-"$4}')

        echo $mid $trim1sum | cat >> $logs/trim1-reads.txt

        fwd_trim2="${trim2}/${mid}-{name}-R1.fq.gz"
        rev_trim2="${trim2}/${mid}-{name}-R2.fq.gz"

        cutadapt -e 0.1 --no-indels -j $1 --discard-untrimmed \
        -g ^file:scratch/rev-adapters.fa \
        -o $rev_trim2 -p $fwd_trim2 $rev_trim1_in $fwd_trim1_in > $logs/trim2-cutadapt.txt

    fi

done

# Apply fastqc to each file with reads after the first trimming step and make multiqc summaries ####
fastqc $(cat $trim1/files.txt | tr '\n' ' ') -o $trim1 -t $1

multiqc -f -o $logs -n trim1-R1.html -ip \
-i 'Forward reads after the first trim' $trim1/*R1_fastqc.zip
multiqc -f -o $logs -n trim1-R2.html -ip \
-i 'Reverse reads after the first trim' $trim1/*R2_fastqc.zip

# Initiate the read count summary for the second and third trimming steps ####
echo trim2 reads | cat >> $logs/trim2-reads.txt
echo trim3 reads | cat >> $logs/trim3-reads.txt

# This last trimming is done to remove any remaining primer sequence in the reverse complement ####
for fwd_trim2_in in $(find $trim2 -name "*rev-R1.fq.gz"); do

    rev_trim2_in=$(echo $fwd_trim2_in | sed 's/R1.fq.gz/R2.fq.gz/')

    #Count sequences after the second trim step
    if [[ -f $fwd_trim2_in && -f $rev_trim2_in ]]; then
        trim2sum=$(echo $(gzip -cd $fwd_trim2_in | wc -l) / 4 | bc)
    else
        trim2sum=0
    fi

    # Trim any remaining primer overhang from reads in files
    if [ $trim2sum -gt 0 ]; then

        echo $fwd_trim2_in | cat >> $trim2/files.txt
        echo $rev_trim2_in | cat >> $trim2/files.txt

        sample=$(basename $fwd_trim2_in | awk -F '-' '{print $1"-"$2"-"$3"-"$4"-"$5}')

        echo $sample $trim2sum | cat >> $logs/trim2-reads.txt

        fwd_trim3="${out}/${sample}-rc-R1.fq.gz"
        rev_trim3="${out}/${sample}-rc-R2.fq.gz"

        if [[ "$(echo $sample | awk '/fun_fwd/ && /fun_rev/' )" == "$sample" ]]; then

            atropos -e 0.1 --aligner insert --insert-match-error-rate 0.2 -T $1 \
            -a $fun_rev_rc -A $fun_fwd_rc \
            -o $fwd_trim3 -p $rev_trim3 -pe1 $fwd_trim2_in -pe2 $rev_trim2_in > $logs/trim3-atropos-fun.txt

        fi

        if [[ "$(echo $sample | awk '/gi_fwd/ && /gi_rev/' )" == "$sample" ]]; then

            atropos -e 0.1 --aligner insert --insert-match-error-rate 0.2 -T $1 \
	    -a $gi_rev_rc -A $gi_fwd_rc \
            -o $fwd_trim3 -p $rev_trim3 -pe1 $fwd_trim2_in -pe2 $rev_trim2_in > $logs/trim3-atropos-gi.txt

        fi
        
        # cutadapt -e 0.1 -j $1 \
        # -a $fun_rev_rc -a $gi_rev_rc -A $fun_fwd_rc -A $gi_fwd_rc \
        # -o $fwd_trim3 -p $rev_trim3 $fwd_trim2_in $rev_trim2_in > $logs/trim3-cutadapt.txt

        if [[ -f $fwd_trim3 && -f $rev_trim3 ]]; then
            trim3sum=$(echo $(gzip -cd $fwd_trim3 | wc -l) / 4 | bc)
        else
            trim3sum=0
        fi

        if [ $trim3sum -gt 0 ]; then
            echo $sample $trim3sum | cat >> $logs/trim3-reads.txt
        fi

    fi

done

# Apply fastqc to each file with reads after the second trimming step and make multiqc summaries ####
fastqc $(cat $trim2/files.txt | tr '\n' ' ') -o $trim2 -t $1

multiqc -f -o $logs -n trim2-R1.html -ip \
-i 'Forward reads after the second trim' $trim2/*R1_fastqc.zip
multiqc -f -o $logs -n trim2-R2.html -ip \
-i 'Reverse reads after the second trim' $trim2/*R2_fastqc.zip

# Apply fastqc to each file with reads after the third trimming step and make multiqc summaries ####
fastqc $(find $out -name "*fq.gz") -o $trim3 -t $1

multiqc -f -o $logs -n trim3-R1.html -ip \
-i 'Forward reads after the third trim' $trim3/*R1_fastqc.zip
multiqc -f -o $logs -n trim3-R2.html -ip \
-i 'Reverse reads after the third trim' $trim3/*R2_fastqc.zip

# Determine how many undetermined reads contain primer sequences ####
cutadapt -e 0.2 -j $1 --discard-untrimmed \
    -g "$fun_fwd;o=${#fun_fwd}" -g "$gi_fwd;o=${#gi_fwd}" -G "$fun_rev;o=${#fun_rev}" -G "$gi_rev;o=${#gi_rev}" \
    -a "$fun_rev_rc;o=${#fun_rev}" -a "$gi_rev_rc;o=${#gi_rev}" -A "$fun_fwd_rc;o=${#fun_fwd}" -A "$gi_fwd_rc;o=${#gi_fwd}" \
    -o $error/trim-R1.fq.gz -p $error/trim-R2.fq.gz $in/undetermined-R1.fq.gz $in/undetermined-R2.fq.gz > $logs/trim-error-cutadapt.txt

echo total trimmed | cat >> $logs/trim-error-reads.txt
total=$(echo $(gzip -cd $in/undetermined-R1.fq.gz | wc -l) / 4 | bc)
trimmed=$(echo $(gzip -cd $error/trim-R1.fq.gz | wc -l) / 4 | bc)
echo $total $trimmed | cat >> $logs/trim-error-reads.txt

# Write fastqc and multiqc reports if any undetermined reads were trimmed ####
if [ $trimmed -gt 0 ]; then
    fastqc $(find $error -name "*fq.gz") -o $error -t $1
    
    multiqc -f -o $logs -n trim-error-R1.html -ip \
    -i 'Undetermined forward reads trimmed with Fun + Gi primers' $error/*R1_fastqc.zip
    
    multiqc -f -o $logs -n trim-error-R2.html -ip \
    -i 'Undetermined reverse reads trimmed with Fun + Gi primers' $error/*R2_fastqc.zip
fi

# Remove the scratch directory after all intermediate logs have been generated ####    
rm -r scratch
