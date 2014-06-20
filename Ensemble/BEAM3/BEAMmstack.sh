#!/bin/bash



cnt=0
fpnum=0
fpcounter=0
tpcounter=0
tp=0


> "$1/finalmresults.txt"

> "$1/finalmresults.match.txt"

echo $cnt > "$1/finalmresults.match.txt"
echo $fpcounter > "$1/finalmfpresults.match.txt"
echo $tfpcounter > "$1/finalmtfpresults.match.txt"
echo $tp > "$1/finalmtpresults.match.txt"

find $1 -type f -name "results*" |
while read file;
do

		echo ""	
		echo "***************************************************"
		echo ""

> "$1/temp.txt"
> "$1/temp2.txt"


		
		awk -F '\t' '{print $4,$3}' $file > "$1/temp.txt"

		sort -g -r "$1/temp.txt" > "$1/temp2.txt"


		check=0;
		ni=0;
		index=0;
		filecnt=0;
		filefpnum=0;
		while read line
		do

			chival=$(cut -d " " -f1 <<<"${line}")
			lbound=17.5
			#echo $chival;
			res=$(echo $chival'>'$lbound | bc -l)
			if [ $res -eq 0 ]; then #33.15 is the limit equal to p-value = 0.05. If greater, p-value is smaller therefore should be accepted
				break
			fi
			echo $line;

			if [[  $line == *" $2"* ]] || [[ $line == *" $3"* ]];
			then
				if [[ "$ni" -eq 1 ]]
					then
						ni=$((ni+1))
						if [[ "$check" -le 0 ]]
						then
						check=$((check+1))
						cnt=$((cnt+1))
						filecnt=$((filecnt+1))
						sed -i '1i\'"$file" "$1/finalmresults.match.txt"
						sed -i '1i\'"$cnt" "$1/finalmresults.match.txt"
						fi
					else
						ni=$((ni+1))
					fi
			else
				fpnum=$((fpnum+1))
				filefpnum=$((filefpnum+1))
			fi

		done < "$1/temp2.txt"
		
		if [ $filefpnum -ge 1 ]; then
			fpcounter=$((fpcounter+1))
			sed -i '1i\'"$file" "$1/finalmfpresults.match.txt"
			sed -i '1i\'"$fpcounter" "$1/finalmfpresults.match.txt"
			sed -i '1i\'"$file" "$1/finalmtfpresults.match.txt"
			sed -i '1i\'"$fpnum" "$1/finalmtfpresults.match.txt"

		else
			tp=$((tp+1))
			sed -i '1i\'"$file" "$1/finalmtpresults.match.txt"
			sed -i '1i\'"$tp" "$1/finalmtpresults.match.txt"
		fi


		rm "$1/temp.txt"
		rm "$1/temp2.txt"

		echo ""	
		echo "$file"
		echo "TP: $filecnt | FP: $filefpnum"
		echo "***************************************************"
		echo ""

		#sed -i '1i\'"Accuracy=$cnt%" "$1/finalresults.txt"
done > "$1/finalmresults.txt"

fcnt=$(sed -n 1p "$1/finalmresults.match.txt")
fpcnt=$(sed -n 1p "$1/finalmfpresults.match.txt")
tpcounter=$(sed -n 1p "$1/finalmtpresults.match.txt")
tfpcounter=$(sed -n 1p "$1/finalmtfpresults.match.txt")
fneg=$((100-fcnt))
gtp=$((fcnt*tpcounter/100))

sed -i '1i\'"FN=$fneg%" "$1/finalmresults.txt"
sed -i '1i\'"FP=$fpcnt%" "$1/finalmresults.txt"
sed -i '1i\'"TP=$fcnt% - $fcnt interactions" "$1/finalmresults.txt"
sed -i '1i\'"Total FP=$tfpcounter interactions" "$1/finalmresults.txt"
sed -i '1i\'"Only TP=$tpcounter% - $gtp% globally ($2 $3)" "$1/finalmresults.txt"

sed -i '1i\'"$fcnt	$fpcnt" "$1/finalmresults.txt"


chmod go+rwx "$1/finalmresults.txt"
chmod go+rwx "$1/finalmtpresults.match.txt"
chmod go+rwx "$1/finalmfpresults.match.txt"
chmod go+rwx "$1/finalmtfpresults.match.txt"
chmod go+rwx "$1/finalmresults.match.txt"