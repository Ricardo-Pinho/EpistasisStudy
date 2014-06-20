#!/bin/bash


oldfile=$1
base=$3
file="$3/$(basename $(dirname ${oldfile}))/$(basename ${oldfile})"
file=${file}.txt

newfile="$2/$(basename $(dirname ${file}))"

temp=$(echo "$file" | sed -e 's/.*1750.1.//g')
ind=$( echo "$temp" | sed -e 's/-3GE-.*//g' )


if [ ! -f "$newfile/results${ind}.txt" ]
then

   mkdir -p "$2/$(basename $(dirname ${file}))"
   chmod go+rxw "$2/$(basename $(dirname ${file}))"


echo " "
echo "*********************************************************"
echo "$ind"
echo "$file"
echo "*********************************************************"
echo " "

( /usr/bin/time -v ./BEAM3 "${file}" -o "results${ind}.txt") 2> "$newfile/time${ind}.txt"


mv "./chi.results${ind}.txt" "$newfile/results${ind}.txt"
mv "./g.results${ind}.txt.dot" "$newfile/g${ind}.dot"
mv "./posterior.results${ind}.txt" "$newfile/posterior${ind}.txt"

chmod 777 "$newfile/results${ind}.txt"
chmod 777 "$newfile/g${ind}.dot"
chmod 777 "$newfile/posterior${ind}.txt"

chmod 777 "$newfile/time${ind}.txt"
#mv "interactions.txt" "interactions${index}.txt"

#chmod 777 "interactions${index}.txt"

rm $file; 

fi