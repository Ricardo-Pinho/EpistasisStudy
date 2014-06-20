#!/bin/bash

oldfile=$1

temp=$(echo "$oldfile" | sed -e 's/.*1750.1.//g')
ind=$(echo "$temp" | sed -e 's/-3GE-.*//g')

mkdir -p "$5/others"
chmod go+rxw "$5/others"

mkdir -p "$5/maineffects/"
chmod go+rxw "$5/maineffects/"



##############################################################################################################################
#BOOST
##############################################################################################################################
> "$5/others/boostMresult${ind}.txt"
base=$2
file="$2/$(basename $(dirname ${oldfile}))/$(basename ${oldfile})"
newfilefolder="$2/$(basename $(dirname ${file}))"

newfile="$newfilefolder/MarginalAssoc${ind}.txt"

#echo "" > "$5/others/boostMresult${ind}.txt"

#while read file;
#do

		#echo ""	>> "$5/others/boostMresult${ind}.txt"
		#echo "***************************************************" >> "$5/others/boostMresult${ind}.txt"
		#echo "$newfile" >> "$5/others/boostMresult${ind}.txt"
		#echo "" >> "$5/others/boostMresult${ind}.txt"

> "$newfilefolder/temp.txt"


		sort -k2 -g -r "$newfile" > "$newfilefolder/temp.txt"



		while read line
		do

			chival=$(cut -d " " -f2 <<<"${line}")
			lbound=17.4
			#echo $line;
			#echo $chival;
			res=$(echo $chival'>'$lbound | bc -l)
			if [ $res -eq 0 ]; then #33.15 is the limit equal to p-value = 0.05. If greater, p-value is smaller therefore should be accepted
				break
			fi

			echo $line

		done < "$newfilefolder/temp.txt" >> "$5/others/boostMresult${ind}.txt"


		rm "$newfilefolder/temp.txt"


		echo "***************************************************" >> "$5/others/boostMresult${ind}.txt"
		#echo "" >> "$5/others/boostMresult${ind}.txt"

		#sed -i '1i\'"Accuracy=$cnt%" "$1/finalresults.txt"
#done < $newfile > "$5/others/boostMresult${ind}.txt"

chmod go+rwx "$5/others/boostMresult${ind}.txt"

##############################################################################################################################
#BEAM 3
##############################################################################################################################
> "$5/others/beam3result${ind}.txt"
base=$3
file="$3/$(basename $(dirname ${oldfile}))/$(basename ${oldfile})"
newfilefolder="$3/$(basename $(dirname ${file}))"

newfile="$newfilefolder/results${ind}.txt"


		#echo ""	>> "$5/others/beam3result${ind}.txt"
		#echo "***************************************************" >> "$5/others/beam3result${ind}.txt"
		#echo "$newfile" >> "$5/others/beam3result${ind}.txt"
		#echo "" >> "$5/others/beam3result${ind}.txt"

> "$newfilefolder/temp.txt"
> "$newfilefolder/temp2.txt"

		awk -F '\t' '{print $4,$3}' $newfile > "$newfilefolder/temp.txt"

		sort -g -r "$newfilefolder/temp.txt" > "$newfilefolder/temp2.txt"

		while read line
		do
			echo $line;
			chival=$(cut -d " " -f1 <<<"${line}")
			lbound=17.5
			#echo $chival;
			res=$(echo $chival'>'$lbound | bc -l)
			if [ $res -eq 0 ]; then #33.15 is the limit equal to p-value = 0.05. If greater, p-value is smaller therefore should be accepted
				break
			fi

		done < "$newfilefolder/temp2.txt" >> "$5/others/beam3result${ind}.txt"

		rm "$newfilefolder/temp2.txt"
		rm "$newfilefolder/temp.txt"


		echo "***************************************************" >> "$5/others/beam3result${ind}.txt"
		#echo "" >> "$5/others/beam3result${ind}.txt"

		#sed -i '1i\'"Accuracy=$cnt%" "$1/finalresults.txt"
#done < $newfile > "$5/others/beam3result${ind}.txt"

chmod go+rwx "$5/others/beam3result${ind}.txt"

##############################################################################################################################
#SNPHarvester
##############################################################################################################################
> "$5/others/snpharvesterresult${ind}.txt"
base=$4
file="$4/$(basename $(dirname ${oldfile}))/$(basename ${oldfile})"
newfilefolder="$4/$(basename $(dirname ${file}))"

newfile="$newfilefolder/results${ind}.txt"


		#echo ""	>> "$5/others/snpharvesterresult${ind}.txt"
		#echo "***************************************************" >> "$5/others/snpharvesterresult${ind}.txt"
		#echo "$newfile" >> "$5/others/snpharvesterresult${ind}.txt"
		#echo "" >> "$5/others/snpharvesterresult${ind}.txt"

> "$newfilefolder/temp.txt"
> "$newfilefolder/temp2.txt"

		awk '/Runnning Time:/{p=0}p;/chi_square/' $newfile > "$newfilefolder/temp.txt"

		while read line
		do
			echo  "${line#*:}"
		done < "$newfilefolder/temp.txt" > "$newfilefolder/temp2.txt"

		sort -r -g "$newfilefolder/temp2.txt" > "$newfilefolder/temp.txt"



		while read line
		do

			if [[ ! $line == *","* ]];
			then
				echo $line;
			fi

done < "$newfilefolder/temp.txt" >> "$5/others/snpharvesterresult${ind}.txt"

rm "$newfilefolder/temp.txt"
rm "$newfilefolder/temp2.txt"


		echo "***************************************************" >> "$5/others/snpharvesterresult${ind}.txt"
		#echo "" >> "$5/others/snpharvesterresult${ind}.txt"

		#sed -i '1i\'"Accuracy=$cnt%" "$1/finalresults.txt"
#done < $newfile > "$5/others/snpharvesterresult${ind}.txt"

chmod go+rwx "$5/others/snpharvesterresult${ind}.txt"

##############################################################################################################################
#Voting
##############################################################################################################################

> "$5/maineffects/maineffectresult${ind}.txt"
echo "***************************************************" >> "$5/maineffects/maineffectresult${ind}.txt"
echo "$ind" >> "$5/maineffects/maineffectresult${ind}.txt"
#echo "" >> "$5/maineffects/maineffectresult${ind}.txt"
#echo "" > "$5/maineffects/maineffectresult${ind}.match.txt"

maineffects="";

while read line
do
	if [[  $line == *"***************************************************"* ]];
	then
		break
	fi
	#echo "$ind"
	n1=$(cut -d " " -f1 <<<"${line}")
	n1r=$((n1+1))
	if  grep -Fq "(rs${n1r}_${n1})"  "$5/others/snpharvesterresult${ind}.txt";
	then
			if [[ ! $maineffects == *"$n1r"* ]];
			then
				maineffects="$maineffects
$n1r";
			fi
	elif  grep -Fq " ${n1r}"  "$5/others/beam3result${ind}.txt";
	then
		harvester=$(head -n 1 "$5/others/snpharvesterresult${ind}.txt")
		if [[ $harvester == *"***************************************************"* ]] && [[ ! $maineffects == *"$n1r"* ]];
		then
				maineffects="$maineffects
$n1r";
		fi
	fi

done < "$5/others/boostMresult${ind}.txt" #>> "$5/maineffects/maineffectresult${ind}.txt"

#echo "A$ind" >> "$5/maineffects/maineffectresult${ind}.txt"

while read line
do
	#echo "X$ind"
	if [[  $line == *"***************************************************"* ]];
	then
		break
	fi
	#echo "Y$ind"
	n1=$(cut -d " " -f2 <<<"${line}")
	n1r=$((n1-1))
	if  grep -Fq "(rs${n1}_${n1r})"  "$5/others/snpharvesterresult${ind}.txt";
	then
			if [[ ! $maineffects == *"$n1"* ]];
			then
				maineffects="$maineffects
$n1r";
			fi
	fi

done < "$5/others/beam3result${ind}.txt" #>> "$5/maineffects/maineffectresult${ind}.txt"

#echo "F$ind" >> "$5/maineffects/maineffectresult${ind}.txt"

echo "$maineffects" >> "$5/maineffects/maineffectresult${ind}.txt"


echo "***************************************************" >> "$5/maineffects/maineffectresult${ind}.txt"
#echo "" >> "$5/maineffects/maineffectresult${ind}.txt"

chmod go+rwx "$5/maineffects/maineffectresult${ind}.txt"