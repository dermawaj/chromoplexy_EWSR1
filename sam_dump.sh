while IFS= read -r row; do
   IFS=" " read VAR1 VAR2 VAR3 VAR4 VAR5 VAR6 VAR7 VAR8 VAR9 EMPTY<<< $row
   echo $VAR1
   file=./dbGAP/"$VAR1"/"$VAR1".sra
   echo $file
   sam-dump $file | samtools view -bS - > ./dbGAP/"$VAR1"/"$VAR1".bam
done < ./SRR_Acc_List.txt