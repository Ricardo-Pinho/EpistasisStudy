#!/bin/bash

oldfile=$1

temp=$(echo "$oldfile" | sed -e 's/.*1750.1.//g')
ind=$(echo "$temp" | sed -e 's/-3GE-.*//g')

mkdir -p "$5/others"
chmod go+rxw "$5/others"

mkdir -p "$5/epistasis"
chmod go+rxw "$5/epistasis"

##############################################################################################################################
#BOOST
##############################################################################################################################
> "$5/others/boostEresult${ind}.txt"
base=$2
file="$2/$(basename $(dirname ${oldfile}))/$(basename ${oldfile})"
newfilefolder="$2/$(basename $(dirname ${file}))"

newfile="$newfilefolder/InteractionRecords${ind}.txt"

echo "" > "$5/others/boostEresult${ind}.txt"

#while read file;
#do

		#echo ""	>> "$5/others/boostEresult${ind}.txt"
		#echo "***************************************************" >> "$5/others/boostEresult${ind}.txt"
		#echo "$newfile" >> "$5/others/boostEresult${ind}.txt"
		#echo "" >> "$5/others/boostEresult${ind}.txt"

> "$newfilefolder/temp.txt"
> "$newfilefolder/temp2.txt"


		awk -F " " '{print $6,$2, $3}' $newfile > "$newfilefolder/temp.txt"

		sort -t " " -g -r "$newfilefolder/temp.txt" > "$newfilefolder/temp2.txt"



		while read line
		do

			chival=$(cut -d " " -f1 <<<"${line}")
			lbound=33.15
			#echo $line
			#echo $chival
			res=$(echo $chival'>'$lbound | bc -l)
			if [ $res -eq 0 ]; then #33.15 is the limit equal to p-value = 0.05. If greater, p-value is smaller therefore should be accepted
				break
			fi

			echo $line

		done < "$newfilefolder/temp2.txt" >> "$5/others/boostEresult${ind}.txt"


		rm "$newfilefolder/temp2.txt"
		rm "$newfilefolder/temp.txt"


		echo "***************************************************" >> "$5/others/boostEresult${ind}.txt"
		#echo "" >> "$5/others/boostEresult${ind}.txt"

		#sed -i '1i\'"Accuracy=$cnt%" "$1/finalresults.txt"
#done < $newfile > "$5/others/boostEresult${ind}.txt"

chmod go+rwx "$5/others/boostEresult${ind}.txt"

##############################################################################################################################
#TEAM
##############################################################################################################################
> "$5/others/teamresult${ind}.txt"
base=$3
file="$3/$(basename $(dirname ${oldfile}))/$(basename ${oldfile})"
newfilefolder="$3/$(basename $(dirname ${file}))"

newfile="$newfilefolder/results${ind}.txt"


		#echo ""	>> "$5/others/teamresult${ind}.txt"
		#echo "***************************************************" >> "$5/others/teamresult${ind}.txt"
		#echo "$newfile" >> "$5/others/teamresult${ind}.txt"
		echo "" > "$5/others/teamresult${ind}.txt"

if  grep -Fxq 'Significant SNP-pair(s) found!'  $newfile;
	then
	
		awk '/================/{p=0}p;/======Significant SNP-pairs=====/{p=1}' $newfile > "$newfilefolder/temp.txt"

		while read line
		do
			echo  "${line#*:};${line%:*}"
		done < "$newfilefolder/temp.txt" > "$newfilefolder/temp2.txt"

		sort -r -g "$newfilefolder/temp2.txt" > "$newfilefolder/temp.txt"

		while read line
		do

			if [[  $line == *"No Significant SNP-pair(s)"* ]];
			then
				break
			fi
			chival=$(echo  "${line%;*}")
			lbound=42.45
			res=$(echo $chival'>'$lbound | bc -l)
			if [ $res -eq 0 ]; then #42.45 is the limit equal to p-value = 0.05. If greater, p-value is smaller therefore should be accepted
			      break
			fi
			echo $line;

		done < "$newfilefolder/temp.txt" >> "$5/others/teamresult${ind}.txt"


		rm "$newfilefolder/temp2.txt"
		rm "$newfilefolder/temp.txt"
fi


		echo "***************************************************" >> "$5/others/teamresult${ind}.txt"
		#echo "" >> "$5/others/teamresult${ind}.txt"

		#sed -i '1i\'"Accuracy=$cnt%" "$1/finalresults.txt"
#done < $newfile > "$5/others/teamresult${ind}.txt"

chmod go+rwx "$5/others/teamresult${ind}.txt"

##############################################################################################################################
#SNPRuler
##############################################################################################################################
> "$5/others/snprulerresult${ind}.txt"
base=$4
file="$4/$(basename $(dirname ${oldfile}))/$(basename ${oldfile})"
newfilefolder="$4/$(basename $(dirname ${file}))"

newfile="$newfilefolder/interactions${ind}.txt"


		#echo ""	>> "$5/others/snprulerresult${ind}.txt"
		#echo "***************************************************" >> "$5/others/snprulerresult${ind}.txt"
		#echo "$newfile" >> "$5/others/snprulerresult${ind}.txt"
		echo "" > "$5/others/snprulerresult${ind}.txt"

sort -k3 -r -g "$newfile" > "$newfilefolder/temp.txt"


while read line
do
	if [[  $line == *"***************************************************"* ]];
	then
		break
	fi
	chival=$(cut -f3 <<<"${line}")
	lbound=42.45
	res=$(echo $chival'>'$lbound | bc -l)
	if [ $res -eq 0 ]; then #42.45 is the limit equal to p-value = 0.05. If greater, p-value is smaller therefore should be accepted
	      break
	fi

	echo $line;

done < "$newfilefolder/temp.txt" >> "$5/others/snprulerresult${ind}.txt"

rm "$newfilefolder/temp.txt"


		echo "***************************************************" >> "$5/others/snprulerresult${ind}.txt"
		#echo "" >> "$5/others/snprulerresult${ind}.txt"

		#sed -i '1i\'"Accuracy=$cnt%" "$1/finalresults.txt"
#done < $newfile > "$5/others/snprulerresult${ind}.txt"

chmod go+rwx "$5/others/snprulerresult${ind}.txt"

##############################################################################################################################
#Voting
##############################################################################################################################
> "$5/epistasis/epistasisresult${ind}.txt"
echo "***************************************************" >> "$5/epistasis/epistasisresult${ind}.txt"
echo "$ind" >> "$5/epistasis/epistasisresult${ind}.txt"
#echo "" >> "$5/epistasis/epistasisresult${ind}.txt"
#echo "" > "$5/epistasis/epistasisresult${ind}.match.txt"

interactions="";

while read line
do
	if [[  $line == *"***************************************************"* ]];
	then
		break
	fi
	#echo "$ind"
	n1=$(cut -d " " -f2 <<<"${line}")
	n2=$(cut -d " " -f3 <<<"${line}")
	n1r=$((n1+1))
	n2r=$((n2+1))
	if  grep -Fq "X${n1r} X${n2r}"  "$5/others/snprulerresult${ind}.txt";
	then
			if [[ ! $interactions == *"$n1r $n2r"* ]];
			then
				interactions="$interactions
$n1r $n2r;$n2r $n1r";
			fi
	elif grep -Fq "X${n2r} X${n1r}"  "$5/others/snprulerresult${ind}.txt";
	then
			if [[ ! $interactions == *"$n1r $n2r"* ]];
			then
				interactions="$interactions
$n1r $n2r;$n2r $n1r";
			fi
	elif  grep -Fq "(${n1},${n2})"  "$5/others/teamresult${ind}.txt";
	then
			if [[ ! $interactions == *"$n1r $n2r"* ]];
			then
				interactions="$interactions
$n1r $n2r;$n2r $n1r";
			fi
	elif  grep -Fq "(${n2},${n1})"  "$5/others/teamresult${ind}.txt";
	then
			if [[ ! $interactions == *"$n1r $n2r"* ]];
			then
				interactions="$interactions
$n1r $n2r;$n2r $n1r";
			fi
	fi

done < "$5/others/boostEresult${ind}.txt" #>> "$5/epistasis/epistasisresult${ind}.txt"

#echo "A$ind" >> "$5/epistasis/epistasisresult${ind}.txt"

while read line
do
	#echo "X$ind"
	if [[  $line == *"***************************************************"* ]];
	then
		break
	fi
	#echo "Y$ind"
	temp="${line#*(}" 
	n1="${temp%,*}"
	temp="${line#*,}"
	n2="${temp%)*}"
	n1r=$((n1+1))
	n2r=$((n2+1))
	if  grep -Fq "X${n1r} X${n2r}"  "$5/others/snprulerresult${ind}.txt";
	then
			if [[ ! $interactions == *"$n1r $n2r"* ]];
			then
				interactions="$interactions
$n1r $n2r;$n2r $n1r";
			fi
	elif grep -Fq "X${n2r} X${n1r}"  "$5/others/snprulerresult${ind}.txt";
	then
			if [[ ! $interactions == *"$n1r $n2r"* ]];
			then
				interactions="$interactions
$n1r $n2r;$n2r $n1r";
			fi
	fi

done < "$5/others/teamresult${ind}.txt" #>> "$5/epistasis/epistasisresult${ind}.txt"

#echo "F$ind" >> "$5/epistasis/epistasisresult${ind}.txt"

#while read line
#do
	#echo "Z$ind"
#	if [[  $line == *"***************************************************"* ]];
#	then
#		break
#	fi
#	temp=$(cut -d " " -f1 <<<"${line}")
#	n1="${temp#*X}"
#	temp=$(cut -d " " -f2 <<<"${line}")
#	n2="${temp#*X}"
	#echo "$n1 $n2"
#	if [[ !  $interactions == *"$n1 $n2"* ]];
#	then
#		interactions="$interactions
#$n1 $n2;$n2 $n1";
#	fi
#
#done < "$5/others/epistasisresult${ind}.txt"
#echo "B$ind" >> "$5/epistasis/epistasisresult${ind}.txt"

echo "$interactions" >> "$5/epistasis/epistasisresult${ind}.txt"


echo "***************************************************" >> "$5/epistasis/epistasisresult${ind}.txt"
#echo "" >> "$5/epistasis/epistasisresult${ind}.txt"

chmod go+rwx "$5/epistasis/epistasisresult${ind}.txt"