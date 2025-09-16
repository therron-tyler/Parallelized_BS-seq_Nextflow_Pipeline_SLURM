# This script assumes your csv and scripts are in the current working directory
# | tr -d '\r' -- put at the end before being written to output script object since is removes  the carriage return character (\r) from the output script. This is necessary because the sbatch command expects Unix line breaks (\n) instead of DOS line breaks (\r\n)

csv=$1
template=$2
script_dir=$(pwd)

template_SN="_SampleName_"
template_FQ="_fastq_"
template_OUT="_output_"

sed 1d $script_dir/$csv | awk '1; END {print ""}' | while read line || [ -n "$line" ]; do
    sample_SN=$(echo $line | cut -d ',' -f 1)
    sample_FQ=$(echo $line | cut -d ',' -f 2)
    sample_OUT=$(echo $line | cut -d ',' -f 3)
    output_script=$script_dir/VDJ_cellranger_$sample_SN.txt
    sed "s:$template_SN:$sample_SN:g" $script_dir/$template | \
    sed "s:$template_FQ:$sample_FQ:g" | 
    sed "s:$template_OUT:$sample_OUT:g" | tr -d '\r' > $output_script
    sbatch $output_script
done
