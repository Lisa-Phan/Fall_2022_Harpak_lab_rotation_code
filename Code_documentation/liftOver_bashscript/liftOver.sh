#!/bin/bash

#export the path  where you install liftover and the chain file
export PATH="$HOME:$PATH"

for i in <list_of_files_to_lift>	
do

#Create a new.BED using chromosome and position columns
#position is duplicated due to SNP having length = 1
awk '{print $2 " " $3 " " $3}' $i > new.BED

#change formatting for column, adding 'chr' in front 
#of all the entries
sed -i -e 's/^/chr/' new.BED

#remove header
sed -i '1d' new.BED

#call liftOver
# liftOver <inputfile> <chainfile> <outputliftfile> <unliftedfile>
liftOver new.BED $HOME/hg19ToHg38.over.chain.gz /scratch/09069/dhp563/Summstats_redo_2/liftOver_fi/conv_$i.BED /scratch/09069/dhp563/Summstats_redo_2/liftOver_fi/unMap_$i.BED

#add a new column header for both out files
sed -i '1i NewChr\tNewPos\tNewPos' /scratch/09069/dhp563/Summstats_redo_2/liftOver_fi/conv_$i.BED
sed -i '1i NewChr\tNewPos\tNewPos' /scratch/09069/dhp563/Summstats_redo_2/liftOver_fi/unMap_$i.BED

#Remove comments (lines with # ) from unMap files
paste -d'\t' $i /scratch/09069/dhp563/Summstats_redo_2/liftOver_fi/conv_$i.BED
sed -i '/#/d' /scratch/09069/dhp563/Summstats_redo_2/liftOver_fi/unMap_$i.BED
done
