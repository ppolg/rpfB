#! /bin/bash -p

bam_dir=$1

outputextension=".count"

shopt -s nullglob
for bam_path in "$bam_dir"/*.bam; do
    bam_file=${bam_path##*/}
    bam_name=${bam_file%.bam}
    sort_name="$bam_name""_sorted.bam"
    echo "Calculating number of reads for $bam_file"
    samtools sort -o "$sort_name" "$bam_file"
    samtools index "$sort_name"
    # htseq-count -f bam -m union -t gene -i ID -a 0 "$bam_file" "$gff_file" > "$bam_name$outputextension"
    echo "saved as "$bam_name$outputextension""
done
echo "All done!"
