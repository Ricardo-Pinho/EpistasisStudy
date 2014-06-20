#!/bin/bash

oldfile=$1

temp=$(echo "$oldfile" | sed -e 's/.*1750.1.//g')
ind=$(echo "$temp" | sed -e 's/-3GE-.*//g')

resultfolder="$3/Ensemble/$(basename $(dirname ${oldfile}))"

mkdir -p "$resultfolder"
chmod go+rxw "$resultfolder"

#echo "$resultfolder"

if [[ ! -f  "$resultfolder/time${ind}.txt" ]] || [[ ! -f  "$resultfolder/epistasis/epistasisresult${ind}.txt" ]] || [[ ! -f  "$resultfolder/maineffects/maineffectresult${ind}.txt" ]];
then

( /usr/bin/time -v bash Ensembleconfig.sh $1 $2 $3 $resultfolder ) 2> "$resultfolder/time${ind}.txt"


chmod go+rxw "$resultfolder/time${ind}.txt"

fi