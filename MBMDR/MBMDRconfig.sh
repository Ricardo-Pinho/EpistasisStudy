#!/bin/bash


oldfile=$1
base=$3
file="$3/$(basename $(dirname ${oldfile}))/$(basename ${oldfile})"
file=${file}.csv
folderPath=$4;
#folderNo=$5;
#folderNo2=$6;

mkdir -p "$2/$(basename $(dirname ${file}))"
chmod go+rxw "$2/$(basename $(dirname ${file}))"
newfile="$2/$(basename $(dirname ${file}))"
chgfile="$(dirname ${file})";

#directory=$( ls -1 $folderPath | sort -n | head -n $folderNo2 | tail -n $folderNo )

if [[ ! -f "$newfile/results${ind}.RData"  ]] #&& [[ $directory == *$(basename $(dirname ${file}))* ]]
then

temp=$(echo "$file" | sed -e 's/.*1750.1.//g')
ind=$( echo "$temp" | sed -e 's/-3GE-.*//g' )
objectfile="$newfile/results${ind}.RData";

( /usr/bin/time -v R --vanilla --file=mbmdrscript.R --args $file $objectfile > "$newfile/results${ind}.txt" ) 2> "$newfile/time${ind}.txt"

chmod 777 "$newfile/results${ind}.txt"
chmod 777 "$newfile/results${ind}.RData"

#mv "interactions.txt" "interactions${index}.txt"

#chmod 777 "interactions${index}.txt"

echo " "
echo "*********************************************************"
echo "$newfile/results${ind}.txt"
echo "*********************************************************"
echo " "

fi