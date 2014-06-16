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

( /usr/bin/time -v java -Xmx7000M -jar SNPHarvester5.jar "${file}" > "$newfile/results${ind}.txt" ) 2> "$newfile/time${ind}.txt"


chmod 777 "$newfile/results${ind}.txt"
chmod 777 "$newfile/time${ind}.txt"

index=$((index+1))
fi