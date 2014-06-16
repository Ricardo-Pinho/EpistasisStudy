#!/bin/bash



oldfile=$1
base=$3
file="$3/$(basename $(dirname ${oldfile}))/$(basename ${oldfile})"
file=${file}.genotype.txt

newfile="$2/$(basename $(dirname ${file}))"

temp=$(echo "$file" | sed -e 's/.*1750.1.//g')
ind=$( echo "$temp" | sed -e 's/-3GE-.*//g' )

if [ ! -f "$newfile/results${ind}.txt" ]
then

	mkdir -p "$2/$(basename $(dirname ${file}))"
	chmod go+rxw "$2/$(basename $(dirname ${file}))"

phenotype=${file:0:${#a}-13};

popnum="$( echo "$file" | sed -e 's#.*caco##; s#.mdr.genotype.txt$##' )"

( /usr/bin/time -v ./team.sh "$file" "$phenotype.phenotype.txt"  "$popnum" 300 100 1 > "$newfile/results${ind}.txt" ) 2> "$newfile/time${ind}.txt"

chmod 777 "$newfile/results${ind}.txt"
chmod 777 "$newfile/time${ind}.txt"

#mv "interactions.txt" "interactions${index}.txt"

#chmod 777 "interactions${index}.txt"

echo " "
echo "*********************************************************"
echo "$file"
echo "$ind"
echo "*********************************************************"
echo " "
fi