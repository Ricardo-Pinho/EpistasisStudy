#!/bin/bash


index=1;
chgfile="test";
file=$1

newfile="$2/$(basename $(dirname ${file}))/$(basename ${file})"

if [ ! -f "${newfile}.genotype.txt" ] || [ ! -f "${newfile}.phenotype.txt" ]
then

if [[ ! $file == *$chgfile* ]];
then
	mkdir -p "$2/$(basename $(dirname ${file}))"
	chgfile="$(dirname ${file})";
	chmod go+rxw "$2/$(basename $(dirname ${file}))"
fi

fchar="";
while read a;
do
#line = ${a:0:${#a}}
linex=${a:4:${#a}};
echo "$linex";
done < $file > ${newfile}.temp.txt

while read a;
do
fchar=${a:0:1};
echo -n "$fchar"
#line = ${a:0:${#a}}
done < $file > ${newfile}.phenotype.txt
index=$((index+1))
echo "$index";

tr -s '[:blank:]' ' ' < ${newfile}.temp.txt > ${newfile}.temp2.txt;

rm ${newfile}.temp.txt;

chmod go+rwx ${newfile}.temp2.txt

./transpose.awk "${newfile}.temp2.txt" > "${newfile}.temp.txt"

rm ${newfile}.temp2.txt;

tr -d '[:blank:]' < ${newfile}.temp.txt > ${newfile}.genotype.txt;
rm ${newfile}.temp.txt;

chmod go+rwx ${newfile}.genotype.txt
chmod go+rwx ${newfile}.phenotype.txt
fi