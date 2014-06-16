#!/bin/bash



cnt=0
fpnum=0
fpcounter=0
tpcounter=0
tp=0

order1="$2 $3"
order2="$3 $2"

> "$1/finalresults.txt"

> "$1/finalresults.match.txt"

echo $cnt > "$1/finalresults.match.txt"
echo $fpcounter > "$1/finalfpresults.match.txt"
echo $tfpcounter > "$1/finaltfpresults.match.txt"
echo $tp > "$1/finaltpresults.match.txt"

find $1 -type f -name "results*.txt" |
while read file;
do

echo ""	
echo "***************************************************"


> "$1/temp.txt"

		sed '1,/ sc1$final/d' $file > "$1/temp.txt"
		#sed -n '/column      snp1 snp2   Estimate Std..Error    z.value     Pr...z.. type/,//p' "$1/temp.txt" > "$1/temp2.txt"
		re='^[0-9]+$'
		#tail -n +2 "$1/temp2.txt"
		check=0;
		n1i=0;
		n2i=0;
		filecnt=0;
		filefpnum=0;
		while read line
		do
			fchar=${line:0:1};
			if [[ $fchar =~ $re ]] && [[ $fchar -ge 2 ]]; then
				echo $line;
				if [[  $line == *$order1* ]] || [[ $line == *$order2* ]];
				then
					if [[ "$check" -le 0 ]]
						then
						check=$((check+1))
						cnt=$((cnt+1))
						filecnt=$((filecnt+1))
						sed -i '1i\'"$file" "$1/finalresults.match.txt"
						sed -i '1i\'"$cnt" "$1/finalresults.match.txt"
						fi
				elif [[  $line == *"$2 <NA>"* ]] || [[  $line == *"<NA> $2"* ]];
				then
						if [[ "$n1i" -eq 0 ]]
						then
							n1i=$((n1i+1))
						fi
				elif [[ $line == *"$3 <NA>"* ]] || [[  $line == *"<NA> $3"* ]];
				then
						if [[ "$n2i" -eq 0 ]]
						then
							n2i=$((n2i+1))
						fi
				else
					fpnum=$((fpnum+1))
					filefpnum=$((filefpnum+1))
				fi
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

	echo "$file"
	echo "TP: $filecnt | FP: $filefpnum"
	echo "***************************************************"
	echo ""

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