#!/bin/bash

MAF001="rs112_111" #267
MAFB001="rs267_266"
MAF005="rs4_3" #239
MAFB005="rs239_238"
MAF01="rs135_134" #230
MAFB01="rs230_229"
MAF03="rs197_196" #266
MAFB03="rs266_265"
MAF05="rs80_79" #229
MAFB05="rs229_228"

echo "" > "$1/tableresults.txt"
echo "" > "$1/tabletresults.txt"
find $1/* -type d | while read -r dir
do
	echo "////////////////////////////////////////////////"
	echo "$dir"

    if [[ $dir == *"0.01"* ]]
	then
	echo "0.01"
	echo "////////////////////////////////////////////////"
    bash SNPHarvesterstack.sh "$dir" "$MAF001" "$MAFB001"	
	elif [[ $dir == *"0.05"* ]]; 
	then
	echo "0.05"
	echo "////////////////////////////////////////////////"
	bash SNPHarvesterstack.sh "$dir" "$MAF005" "$MAFB005"
	elif [[ $dir == *"0.1"* ]]; 
	then
	echo "0.1"
	echo "////////////////////////////////////////////////"
	bash SNPHarvesterstack.sh "$dir" "$MAF01" "$MAFB01"
	elif [[ $dir == *"0.3"* ]]; 
	then
	echo "0.3"
	echo "////////////////////////////////////////////////"
	bash SNPHarvesterstack.sh "$dir" "$MAF03" "$MAFB03"
	else
	echo "0.5"
	echo "////////////////////////////////////////////////"
	bash SNPHarvesterstack.sh "$dir" "$MAF05" "$MAFB05"
	fi
	bash SNPHarvestertstack.sh "$dir"
	toverall=$(sed -n 1p "$dir/finaltresults.txt")
	dirname=$(basename ${dir})
	sed -i '1i\'"$dirname		$toverall" "$1/tabletresults.txt"
	chmod go+rwx "$1/tabletresults.txt"

	overall=$(sed -n 1p "$dir/finalresults.txt")
	sed -i '1i\'"$dirname		$overall" "$1/tableresults.txt"
	chmod go+rwx "$1/tableresults.txt"
done