#!/bin/bash


oldfile=$1
base=$3
file="$3/$(basename $(dirname ${oldfile}))/$(basename ${oldfile})"
file=${file}.txt

newfile="$2/$(basename $(dirname ${file}))"

temp=$(echo "$file" | sed -e 's/.*1750.1.//g')
ind=$(echo "$temp" | sed -e 's/-3GE-.*//g')


if [ ! -f "$newfile/InteractionRecords${ind}.txt" ]
then

   mkdir -p "$2/$(basename $(dirname ${file}))"
   chmod go+rxw "$2/$(basename $(dirname ${file}))"

> filenamelist.txt

if [[ "${file}" != *filename* ]];
then
   echo "${file}" >> filenamelist.txt

   ( /usr/bin/time -v ./BOOST ) 2> "$newfile/time${ind}.txt"

   chmod 777 InteractionRecords.txt

   mv "InteractionRecords.txt"  "$newfile/InteractionRecords${ind}.txt"

   chmod 777 MarginalAssoc.txt

   mv "MarginalAssoc.txt"  "$newfile/MarginalAssoc${ind}.txt"

   rm filenamelist.txt

   chmod 777 "$newfile/time${ind}.txt"


   index=$((index+1))
	echo " "
	echo "*********************************************************"
	echo "$ind"
	echo "*********************************************************"
	echo " "
fi

fi

#./BOOST

#chmod 777 InteractionRecords.txt

#chmod 777 MarginalAssoc.txt
