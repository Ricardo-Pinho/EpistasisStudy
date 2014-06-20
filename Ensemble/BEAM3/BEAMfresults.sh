#!/bin/bash

MAF001="112" #267
MAFB001="267"
MAF005="4" #239
MAFB005="239"
MAF01="135" #230
MAFB01="230"
MAF03="197" #266
MAFB03="266"
MAF05="80" #229
MAFB05="229"


echo "" > "$1/tableresults.txt"
echo "" > "$1/tabletresults.txt"
echo "" > "$1/tablemresults.txt"
echo "" > "$1/tablesumresults.txt"
find $1/* -type d | while read -r dir
do
	if [[ $dir == *"ME+I"* ]] || [[ $dir == *"ME"* ]]
	then

	echo "////////////////////////////////////////////////"
	echo "$dir"

    if [[ $dir == *"0.01"* ]]
	then
	echo "0.01"
	echo "////////////////////////////////////////////////"
    bash BEAMmstack.sh "$dir" "$MAF001" "$MAFB001"	
	elif [[ $dir == *"0.05"* ]]; 
	then
	echo "0.05"
	echo "////////////////////////////////////////////////"
	bash BEAMmstack.sh "$dir" "$MAF005" "$MAFB005"
	elif [[ $dir == *"0.1"* ]]; 
	then
	echo "0.1"
	echo "////////////////////////////////////////////////"
	bash BEAMmstack.sh "$dir" "$MAF01" "$MAFB01"
	elif [[ $dir == *"0.3"* ]]; 
	then
	echo "0.3"
	echo "////////////////////////////////////////////////"
	bash BEAMmstack.sh "$dir" "$MAF03" "$MAFB03"
	else
	echo "0.5"
	echo "////////////////////////////////////////////////"
	bash BEAMmstack.sh "$dir" "$MAF05" "$MAFB05"
	fi
	if [[ ! $dir == *"ME+I"* ]]
	then
		dirname=$(basename ${dir})
		overall=$(sed -n 1p "$dir/finalmresults.txt")
		sed -i '1i\'"$dirname		$overall" "$1/tablemresults.txt"
		chmod go+rwx "$1/tablemresults.txt"
	fi
fi

if [[ $dir == *"ME+I"* ]] || [[ $dir == *"I"* ]]
then	
	echo "////////////////////////////////////////////////"
	echo "$dir"

    if [[ $dir == *"0.01"* ]]
	then
	echo "0.01"
	echo "////////////////////////////////////////////////"
    bash BEAMstack.sh "$dir" "$MAF001" "$MAFB001"	
	elif [[ $dir == *"0.05"* ]]; 
	then
	echo "0.05"
	echo "////////////////////////////////////////////////"
	bash BEAMstack.sh "$dir" "$MAF005" "$MAFB005"
	elif [[ $dir == *"0.1"* ]]; 
	then
	echo "0.1"
	echo "////////////////////////////////////////////////"
	bash BEAMstack.sh "$dir" "$MAF01" "$MAFB01"
	elif [[ $dir == *"0.3"* ]]; 
	then
	echo "0.3"
	echo "////////////////////////////////////////////////"
	bash BEAMstack.sh "$dir" "$MAF03" "$MAFB03"
	else
	echo "0.5"
	echo "////////////////////////////////////////////////"
	bash BEAMstack.sh "$dir" "$MAF05" "$MAFB05"
	fi
	if [[ ! $dir == *"ME+I"* ]]
	then
		dirname=$(basename ${dir})
		overall=$(sed -n 1p "$dir/finalresults.txt")
		sed -i '1i\'"$dirname		$overall" "$1/tableresults.txt"
		chmod go+rwx "$1/tableresults.txt"
	fi
fi
if [[ $dir == *"ME+I"* ]]
then
	bash BEAMsumstack.sh "$dir" "$dir/finalmresults.txt" "$dir/finalresults.txt"
	dirname=$(basename ${dir})
	overall=$(sed -n 1p "$dir/finalsumresults.txt")
	sed -i '1i\'"$dirname		$overall" "$1/tablesumresults.txt"
	chmod go+rwx "$1/tablesumresults.txt"
fi
	bash BEAMtstack.sh "$dir"
	toverall=$(sed -n 1p "$dir/finaltresults.txt")
	dirname=$(basename ${dir})
	sed -i '1i\'"$dirname		$toverall" "$1/tabletresults.txt"
	chmod go+rwx "$1/tabletresults.txt"
done