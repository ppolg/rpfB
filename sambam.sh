#! /bin/bash -p

sam_dir=$1

shopt -s nullglob
for sam_path in "$sam_dir"/*.sam; do
    sam_file=${sam_path##*/}
    sam_name=${sam_file%.sam}
    bam_file="$sam_name.bam"
    echo "making $sam_name into a bam."
    samtools view -bo $bam_file $sam_file
done
echo "All done! :)"