#!/bin/bash

MAF001="112" #267
MAFB001="267" #267

MAF005="4" #239
MAFB005="239" #239

MAF01="135" #230
MAFB01="230" #230

MAF03="197" #266
MAFB03="266" #266

MAF05="80" #229
MAFB05="229" #229


echo "" > "$1/tableresults.txt"
echo "" > "$1/tabletresults.txt"
find $1/* -type d | while read dir
do
if [[ ! $dir == *"epistasis"* ]] && [[ ! $dir == *"maineffects"* ]] && [[ ! $dir == *"others"* ]];
then


		echo "////////////////////////////////////////////////"
		echo "$dir"

	    if [[ $dir == *"0.01"* ]]
		then
		echo "0.01"
		echo "////////////////////////////////////////////////"
	    bash Ensemblestack.sh "$dir" "$MAF001" "$MAFB001"	
		elif [[ $dir == *"0.05"* ]]; 
		then
		echo "0.05"
		echo "////////////////////////////////////////////////"
		bash Ensemblestack.sh "$dir" "$MAF005" "$MAFB005"
		elif [[ $dir == *"0.1"* ]]; 
		then
		echo "0.1"
		echo "////////////////////////////////////////////////"
		bash Ensemblestack.sh "$dir" "$MAF01" "$MAFB01"
		elif [[ $dir == *"0.3"* ]]; 
		then
		echo "0.3"
		echo "////////////////////////////////////////////////"
		bash Ensemblestack.sh "$dir" "$MAF03" "$MAFB03"
		else
		echo "0.5"
		echo "////////////////////////////////////////////////"
		bash Ensemblestack.sh "$dir" "$MAF05" "$MAFB05"
		fi
		dirname=$(basename ${dir})
		overall=$(sed -n 1p "$dir/finalresults.txt")
		sed -i '1i\'"$dirname		$overall" "$1/tableresults.txt"
		chmod go+rwx "$1/tableresults.txt"


		bash Ensembletstack.sh "$dir"
		toverall=$(sed -n 1p "$dir/finaltresults.txt")
		dirname=$(basename ${dir})
		sed -i '1i\'"$dirname		$toverall" "$1/tabletresults.txt"
		chmod go+rwx "$1/tabletresults.txt"
fi
done