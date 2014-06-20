#!/bin/bash



cnt=0
fpnum=0
fpcounter=0
tpcounter=0
tp=0

order1="$2 $3"

> "$1/finalresults.txt"

> "$1/finalresults.match.txt"

echo $cnt > "$1/finalresults.match.txt"
echo $fpcounter > "$1/finalfpresults.match.txt"
echo $tfpcounter > "$1/finaltfpresults.match.txt"
echo $tp > "$1/finaltpresults.match.txt"

find "$1/epistasis" -type f -name "epistasisresult*" |
while read file;
do

temp=$(echo "$file" | sed -e 's/.*epistasisresult//g')
ind=$(echo "$temp" | sed -e 's/.txt.*//g')

#resultfolder="$1/Ensemble/$(basename $(dirname ${oldfile}))"

		echo ""	
		echo "***************************************************"
		echo "$ind"
		echo "---------------------------------------------------"


		check=0;
		ni=0;
		index=0;
		filecnt=0;
		filefpnum=0;

		if [[ $1 == *"I"* ]]; then
			while read line
			do
				if [[  $line == *"***************************************************"* ]];
				then
					break
				fi
				echo $line;

				if [[  $line == *$order1* ]]
				then
					if [[ "$check" -le 0 ]]
						then
						check=$((check+1))
						cnt=$((cnt+1))
						filecnt=$((filecnt+1))
						sed -i '1i\'"$file" "$1/finalresults.match.txt"
						sed -i '1i\'"$cnt" "$1/finalresults.match.txt"
						fi
				else
					fpnum=$((fpnum+1))
					filefpnum=$((filefpnum+1))
				fi

			done < <(tail -n "+4" $file)
		fi

		if [[ $1 == *"ME"* ]]; then
			while read line
			do
				if [[  $line == *"***************************************************"* ]];
				then
					break
				fi
				echo $line;

				if [[  $line == *$2* ]] || [[ $line == *$3* ]];
				then
						if [[ "$ni" -eq 1 ]]
						then
							ni=$((ni+1))
							if [[ "$check" -le 0 ]]
							then
							check=$((check+1))
							cnt=$((cnt+1))
							filecnt=$((filecnt+1))
							sed -i '1i\'"$file" "$1/finalresults.match.txt"
							sed -i '1i\'"$cnt" "$1/finalresults.match.txt"
							fi
						else
							ni=$((ni+1))
						fi
				else
					fpnum=$((fpnum+1))
					filefpnum=$((filefpnum+1))
				fi

			done < <(tail -n "+4" "$1/maineffects/maineffectresult${ind}.txt")
		fi

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

		echo "TP: $filecnt | FP: $filefpnum"
		echo "***************************************************"
		echo ""

		#sed -i '1i\'"Accuracy=$cnt%" "$1/finalresults.txt"
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