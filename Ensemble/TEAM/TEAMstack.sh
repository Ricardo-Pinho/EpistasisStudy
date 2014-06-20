#!/bin/bash



cnt=0
fpnum=0
fpcounter=0
tpcounter=0
tp=0

order1="$2,$3"
order2="$3,$2"

> "$1/finalresults.txt"

> "$1/finalresults.match.txt"

echo $cnt > "$1/finalresults.match.txt"
echo $fpcounter > "$1/finalfpresults.match.txt"
echo $tfpcounter > "$1/finaltfpresults.match.txt"
echo $tp > "$1/finaltpresults.match.txt"

find $1 -type f -name "results*" |
while read file;
do

	if  ! grep -Fxq 'Significant SNP-pair(s) found!'  $file;
	then
		echo ""	
		echo "***************************************************"
		echo "$file"
		echo "No Significant SNP-pair(s) found!"
		echo "***************************************************"
		echo ""
	else

		echo ""	
		echo "***************************************************"
		echo "$file"
	
		awk '/================/{p=0}p;/======Significant SNP-pairs=====/{p=1}' $file > "$1/temp.txt"

		index=0;
		filecnt=0;
		filefpnum=0;
		while read line
		do
			echo  "${line#*:};${line%:*}"
		done < "$1/temp.txt" > "$1/temp2.txt"

		sort -r -g "$1/temp2.txt" > "$1/temp.txt"

		index=0;
		while read line
		do

			if [[  $line == *"No Significant SNP-pair(s)"* ]];
			then
				filecnt=0;
				break
			fi
			chival=$(echo  "${line%;*}")
			lbound=42.45
			res=$(echo $chival'>'$lbound | bc -l)
			if [ $res -eq 0 ]; then #42.45 is the limit equal to p-value = 0.05. If greater, p-value is smaller therefore should be accepted
			      break
			fi
			echo $line;

			if [[  $line == *$order1* ]] || [[ $line == *$order2* ]];
			then
			if [[ $index -eq 0 ]]
				then
					index=$((index+1))
					cnt=$((cnt+1))
					filecnt=$((filecnt+1))
					sed -i '1i\'"$file" "$1/finalresults.match.txt"
					sed -i '1i\'"$cnt" "$1/finalresults.match.txt"
				fi
			else
				fpnum=$((fpnum+1))
				filefpnum=$((filefpnum+1))
			fi

		done < "$1/temp.txt" 

		if [ $filefpnum -ge 1 ]; then
			fpcounter=$((fpcounter+1))
			sed -i '1i\'"$file" "$1/finalfpresults.match.txt"
			sed -i '1i\'"$fpcounter" "$1/finalfpresults.match.txt"
			sed -i '1i\'"$file" "$1/finaltfpresults.match.txt"
			sed -i '1i\'"$fpnum" "$1/finaltfpresults.match.txt"

		else
			tp=$((tp+1))
			sed -i '1i\'"$file" "$1/finaltpresults.match.txt"
			sed -i '1i\'"$tp" "$1/finaltpresults.match.txt"
		fi

		rm "$1/temp.txt"
		rm "$1/temp2.txt"

		echo "TP: $filecnt | FP: $filefpnum"
		echo "***************************************************"
		echo ""
		#sed -i '1i\'"Accuracy=$cnt%" "$1/finalresults.txt"
	fi

done > "$1/finalresults.txt"


fcnt=$(sed -n 1p "$1/finalresults.match.txt")
fpcnt=$(sed -n 1p "$1/finalfpresults.match.txt")
tpcounter=$(sed -n 1p "$1/finaltpresults.match.txt")
tfpcounter=$(sed -n 1p "$1/finaltfpresults.match.txt")
fneg=$((100-fcnt))
gtp=$((fcnt*tpcounter/100))

sed -i '1i\'"FN=$fneg%" "$1/finalresults.txt"
sed -i '1i\'"FP=$fpcnt%" "$1/finalresults.txt"
sed -i '1i\'"TP=$fcnt% - $fcnt interactions" "$1/finalresults.txt"
sed -i '1i\'"Total FP=$tfpcounter interactions" "$1/finalresults.txt"
sed -i '1i\'"Only TP=$tpcounter% - $gtp% globally ($2 $3)" "$1/finalresults.txt"

sed -i '1i\'"$fcnt	$fpcnt" "$1/finalresults.txt"


chmod go+rwx "$1/finalresults.txt"
chmod go+rwx "$1/finaltpresults.match.txt"
chmod go+rwx "$1/finalfpresults.match.txt"
chmod go+rwx "$1/finaltfpresults.match.txt"
chmod go+rwx "$1/finalresults.match.txt"