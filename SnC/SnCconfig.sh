#!/bin/bash


oldfile=$1
base=$3
file="$3/$(basename $(dirname ${oldfile}))/$(basename ${oldfile})"
file=${file}.genotype.csv
newfile="$2/$(basename $(dirname ${file}))"

temp=$(echo "$file" | sed -e 's/.*1750.1.//g')
ind=$(echo "$temp" | sed -e 's/-3GE-.*//g')


if [[ ! -f "$newfile/results${ind}.txt"  ]] #|| [[ ! -f "$newfile/results${ind}.RData" ]]
then

   mkdir -p "$2/$(basename $(dirname ${file}))"
   chmod go+rxw "$2/$(basename $(dirname ${file}))"

temp2=$( echo "$file" | sed -e 's/.genotype.csv.*//g' )
phenotype=$( echo "$temp2.phenotype.csv" )

objectfile="$newfile/results${ind}.RData";

( /usr/bin/time -v R --vanilla --file=SnCscript.R --args $file $phenotype $objectfile  > "$newfile/results${ind}.txt" ) 2> "$newfile/time${ind}.txt"

chmod 777 "$newfile/results${ind}.txt"

chmod 777 "$newfile/results${ind}.RData"

chmod 777 "$newfile/time${ind}.txt"

#mv "interactions.txt" "interactions${index}.txt"

#chmod 777 "interactions${index}.txt"

echo " "
echo "*********************************************************"
echo "$file"
echo "$phenotype"
echo "*********************************************************"
echo " "

fi