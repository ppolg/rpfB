#! /bin/bash -p

dir=$1
index=$2
echo "Hi! I am mapping your paired fastq files to your $index index! :)"

shopt -s nullglob
for k in "$dir"/*R1_001.fastq; do
    m=${k/R1_001.fastq/R2_001.fastq}
    kshort=${k##*/}
    mshort=${m##*/}
    outfile=${kshort/_R*}
    bamname="$outfile"".bam"
    samname="$outfile"".sam"
    echo "mapping $kshort and $mshort"
    bowtie2 -x "$index" -1 "$k" -2 "$m" > "$samname"
    echo "saved as $samname"
done
echo "All done! :)"