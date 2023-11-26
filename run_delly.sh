#!/bin/bash

conda activate delly
cd /path_to/dbGAP
#!/bin/bash

# Initialize an empty variable to store the concatenated string
concatenated_str=""

touch samples_control.tsv
# Process each line in the "tumor_normal.txt" file
while IFS= read -r line
do
  # Extract the third column
  str=$(echo "$line" | awk '{print $3}')
  
  # Construct the desired file path
  file_path="./dbGAP/$str/$str.bam"
  
  # Append the file path to the concatenated string, separated by a space
  concatenated_str="$concatenated_str $file_path"

  control_case=$(samtools view -H $file_path | grep -m 1 -Eo '@RG\s+.*\s+SM:([^[:space:]]+)' | grep -o 'SM:[^[:space:]]*'| cut -d':' -f2 | cut -d$'\t' -f1)
  echo -e "$control_case\tcontrol" >> samples_control.tsv
done < tumor_normal.txt

# Remove the leading space (if present) and save the result in a variable
concatenated_str="${concatenated_str# }"


# Print the concatenated string
echo "Concatenated string: $concatenated_str"
# while IFS= read -r row; do
#    IFS=" " read VAR1 tumor normal VAR3 VAR4 VAR5 VAR6 VAR7 VAR8 VAR9 EMPTY<<< $row
#    cd /data/vanderbilt/dermawaj/dbGAP
#    echo $tumor
#    echo $normal
#    tumor_file=./dbGAP/"$tumor"/"$tumor".bam
#    normal_file=./dbGAP/"$normal"/"$normal".bam
#    # echo -e "$tumor\ttumor" > samples.tsv
#    echo $tumor_file
#    echo $normal_file
#    samtools index $tumor_file
#    samtools index $normal_file


# done < ./tumor_normal.txt

# Path: path_to/run_delly.sh

while IFS= read -r row; do
  IFS=" " read VAR1 tumor normal VAR3 VAR4 VAR5 VAR6 VAR7 VAR8 VAR9 EMPTY<<< $row
  cd /data/vanderbilt/dermawaj/dbGAP
  echo $tumor
  echo $normal
  tumor_file=./dbGAP/"$tumor"/"$tumor".bam
  normal_file=./dbGAP/"$normal"/"$normal".bam
  # echo -e "$tumor\ttumor" > samples.tsv
  echo $tumor_file
  echo $normal_file
#  samtools index $tumor_file
#  samtools index $normal_file

  tumor_case=$(samtools view -H $tumor_file | grep -m 1 -Eo '@RG\s+.*\s+SM:([^[:space:]]+)' | grep -o 'SM:[^[:space:]]*'| cut -d':' -f2 | cut -d$'\t' -f1)
  control_case=$(samtools view -H $normal_file | grep -m 1 -Eo '@RG\s+.*\s+SM:([^[:space:]]+)' | grep -o 'SM:[^[:space:]]*'| cut -d':' -f2 | cut -d$'\t' -f1)
  echo -e "$tumor_case\ttumor" > samples.tsv
  echo -e "$control_case\tcontrol" >> samples.tsv
  # echo $tumor_case $control_case
  #pause 1 min 
  # sleep 60
  # delly call -x ./hg19.excl.tsv -o ./delly/"$tumor".bcf -g ./Homo_sapiens_assembly19.fasta $tumor_file $normal_file
  delly call -x ./hg19.excl.tsv -o ./delly/"$tumor".bcf -g /data/vanderbilt/dermawaj/dbGAP/Homo_sapiens_assembly19.fasta $tumor_file $normal_file
  
  delly filter -f somatic -o ./delly/"$tumor".pre.bcf -s samples.tsv ./delly/"$tumor".bcf 
  

  delly call -g /data/vanderbilt/dermawaj/dbGAP/Homo_sapiens_assembly19.fasta -v ./delly/"$tumor".pre.bcf -o ./delly/"$tumor".geno.bcf -x ./hg19.excl.tsv $tumor_file $concatenated_str

  echo -e "$tumor_case\ttumor" > samples_w_control.tsv
  cat samples_control.tsv >> samples_w_control.tsv

  delly filter -f somatic -o ./delly/"$tumor".somatic.bcf -s samples_w_control.tsv ./delly/"$tumor".geno.bcf
  bcftools view --output-type v -o ./delly/"$tumor".pre.vcf  ./delly/"$tumor".pre.bcf 
  bcftools view --output-type v -o ./delly/"$tumor".vcf  ./delly/"$tumor".bcf
  bcftools view --output-type v -o ./delly/"$tumor".somatic.vcf  ./delly/"$tumor".somatic.bcf


done < ./tumor_normal.txt

