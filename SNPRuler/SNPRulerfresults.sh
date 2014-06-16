#!/bin/bash

MAF001="X112" #267
MAFB001="X267"
MAF005="X4" #239
MAFB005="X239"
MAF01="X135" #230
MAFB01="X230"
MAF03="X197" #266
MAFB03="X266"
MAF05="X80" #229
MAFB05="X229"

echo "" > "$1/tableresults.txt"
echo "" > "$1/tabletresults.txt"
find $1/* -type d | while read -r dir
do

if [[ ! $dir == *"ME+I"* ]] && [[ $dir == *"I"* ]]
then
	echo "////////////////////////////////////////////////"
	echo "$dir"
    if [[ $dir == *"0.01"* ]]
	then
	echo "0.01"
	echo "////////////////////////////////////////////////"
    bash SNPRulerstack.sh "$dir" "$MAF001" "$MAFB001"	
	elif [[ $dir == *"0.05"* ]]; 
	then
	echo "0.05"
	echo "////////////////////////////////////////////////"
	bash SNPRulerstack.sh "$dir" "$MAF005" "$MAFB005"
	elif [[ $dir == *"0.1"* ]]; 
	then
	echo "0.1"
	echo "////////////////////////////////////////////////"
	bash SNPRulerstack.sh "$dir" "$MAF01" "$MAFB01"
	elif [[ $dir == *"0.3"* ]]; 
	then
	echo "0.3"
	echo "////////////////////////////////////////////////"
	bash SNPRulerstack.sh "$dir" "$MAF03" "$MAFB03"
	else
	echo "0.5"
	echo "////////////////////////////////////////////////"
	bash SNPRulerstack.sh "$dir" "$MAF05" "$MAFB05"
	fi
	bash SNPRulertstack.sh "$dir"
	toverall=$(sed -n 1p "$dir/finaltresults.txt")
	dirname=$(basename ${dir})
	sed -i '1i\'"$dirname		$toverall" "$1/tabletresults.txt"
	
	overall=$(sed -n 1p "$dir/finalresults.txt")
	dirname=$(basename ${dir})
	sed -i '1i\'"$dirname		$overall" "$1/tableresults.txt"
fi
chmod go+rwx "$1/tableresults.txt"
chmod go+rwx "$1/tabletresults.txt"
done