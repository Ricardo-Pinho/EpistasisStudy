#!/bin/bash


lock1=0;
lock2=0;
dataset1=0;
dataset2=0;
TP1=0
TP2=0
FP1=0
FP2=0
#IFS1=1;
#IFS2=1;
lineA="test"
lineB="test"

> "$1/finalsumresults.txt"
cnt=0
fpnum=0
fpcounter=0
tpcounter=0
tp=0
echo $cnt > "$1/finalsumresults.match.txt"
echo $fpcounter > "$1/finalsumfpresults.match.txt"
echo $tfpcounter > "$1/finalsumtfpresults.match.txt"
echo $tp > "$1/finalsumtpresults.match.txt"

while true; do

if [[ $lock1 -eq 0 ]] && [[ $lock2 -eq 0 ]]; then
	#IFS= read -r lineA <&2
	#IFS= read -r lineB <&3
	#echo $lineA;
	#echo $lineB;
	if ! IFS= read -r lineA <&2 && ! IFS= read -r lineB <&3; then
		echo "wtf"
		break
	fi
else
	if [[ $lock1 -eq 0 ]]; then
		IFS= read -r lineA <&2
	fi
	if [[ $lock2 -eq 0 ]]; then
		IFS= read -r lineB <&3
	fi
fi

  #echo "$lineA"
  #echo "$lineB"
  if [[  $lineA == *"MarginalAssoc"* ]]; then
  		temp=$(echo  "${lineA#*MarginalAssoc}")
  		dataset1=$(echo  "${temp%.*}")
		read -r lineA <&2;
		temp=$(echo  "${lineA#*TP: }")
		TP1=$(echo  "${temp% |*}")
		FP1=$(echo  "${lineA#*FP: }")
		#echo "**$dataset1"
		#echo "TP: $TP1 | FP: $FP1"
		lock1=1;
  		if [[ $lock2 -eq 1 ]]; then
	  		if [[ ! $dataset1 == $dataset2 ]]; then
	  			echo "$dataset1 - $dataset2"
  				echo "error mismatch"
  			fi
  			check=0;
			ni=0;
			index=0;
			filecnt=0;
			filefpnum=0;
  			#echo "$dataset"
  			TP=0;
  			FP=$((FP1+FP2));
  			if [[ $TP1 -eq 1 ]] || [[ $TP2 -eq 1 ]]; then
  				TP=1;
  				cnt=$((cnt+1))
				filecnt=$((filecnt+1))
				sed -i '1i\'"$file" "$1/finalsumresults.match.txt"
				sed -i '1i\'"$cnt" "$1/finalsumresults.match.txt"
  			fi 

  			if [[ $FP -ge 1 ]]; then
				fpnum=$((fpnum+$FP))
				filefpnum=$((filefpnum+1))
  			fi
  			if [ $filefpnum -ge 1 ]; then
			fpcounter=$((fpcounter+1))
			sed -i '1i\'"$file" "$1/finalsumfpresults.match.txt"
			sed -i '1i\'"$fpcounter" "$1/finalsumfpresults.match.txt"
			sed -i '1i\'"$file" "$1/finalsumtfpresults.match.txt"
			sed -i '1i\'"$fpnum" "$1/finalsumtfpresults.match.txt"
			else
				tp=$((tp+1))
				sed -i '1i\'"$file" "$1/finalsumtpresults.match.txt"
				sed -i '1i\'"$tp" "$1/finalsumtpresults.match.txt"
			fi

  			echo ""	
			echo "***************************************************"
			echo "$dataset1"
			echo "Main Effect:"
			echo "TP: $TP1 | FP: $FP1"
			echo "Epistasis:"
			echo "TP: $TP2 | FP: $FP2"
			echo "TP: $TP | FP: $FP"
			echo "***************************************************"
			echo ""
			lock2=0;
			lock1=0;
  		fi
  fi

  if [[  $lineB == *"InteractionRecords"* ]]; then
  		temp=$(echo  "${lineB#*InteractionRecords}")
  		dataset2=$(echo  "${temp%.*}")
		read -r lineB <&3;
		temp=$(echo  "${lineB#*TP: }")
		TP2=$(echo  "${temp% |*}")
		FP2=$(echo  "${lineB#*FP: }")
		#echo "--$dataset2"
		#echo "TP: $TP2 | FP: $FP2"
		lock2=1;
		if [[ $lock1 -eq 1 ]]; then
	  		if [[ ! $dataset1 == $dataset2 ]]; then
				echo "$dataset1 - $dataset2"
				echo "error mismatch"
			fi
			check=0;
			ni=0;
			index=0;
			filecnt=0;
			filefpnum=0;
			#echo "$dataset"
			TP=0;
			FP=$((FP1+FP2));
			if [[ $TP1 -eq 1 ]] || [[ $TP2 -eq 1 ]]; then
				TP=1;
				cnt=$((cnt+1))
				filecnt=$((filecnt+1))
				sed -i '1i\'"$file" "$1/finalsumresults.match.txt"
				sed -i '1i\'"$cnt" "$1/finalsumresults.match.txt"
  			fi 

  			if [[ $FP -ge 1 ]]; then
				fpnum=$((fpnum+$FP))
				filefpnum=$((filefpnum+1))
  			fi
  			if [ $filefpnum -ge 1 ]; then
				fpcounter=$((fpcounter+1))
				sed -i '1i\'"$file" "$1/finalsumfpresults.match.txt"
				sed -i '1i\'"$fpcounter" "$1/finalsumfpresults.match.txt"
				sed -i '1i\'"$file" "$1/finalsumtfpresults.match.txt"
				sed -i '1i\'"$fpnum" "$1/finalsumtfpresults.match.txt"
			else
				tp=$((tp+1))
				sed -i '1i\'"$file" "$1/finalsumtpresults.match.txt"
				sed -i '1i\'"$tp" "$1/finalsumtpresults.match.txt"
			fi
			
			echo ""	
			echo "***************************************************"
			echo "$dataset1"
			echo "Main Effect:"
			echo "TP: $TP1 | FP: $FP1"
			echo "Epistasis:"
			echo "TP: $TP2 | FP: $FP2"
			echo "Total:"
			echo "TP: $TP | FP: $FP"
			echo "***************************************************"
			echo ""
			lock1=0;
			lock2=0;
  		fi
  fi
done 2<$2 3<$3 > "$1/finalsumresults.txt"

fcnt=$(sed -n 1p "$1/finalsumresults.match.txt")
fpcnt=$(sed -n 1p "$1/finalsumfpresults.match.txt")
tpcounter=$(sed -n 1p "$1/finalsumtpresults.match.txt")
tfpcounter=$(sed -n 1p "$1/finalsumtfpresults.match.txt")
fneg=$((100-fcnt))
gtp=$((fcnt*tpcounter/100))

sed -i '1i\'"FN=$fneg%" "$1/finalsumresults.txt"
sed -i '1i\'"FP=$fpcnt%" "$1/finalsumresults.txt"
sed -i '1i\'"TP=$fcnt% - $fcnt interactions" "$1/finalsumresults.txt"
sed -i '1i\'"Total FP=$tfpcounter interactions" "$1/finalsumresults.txt"
sed -i '1i\'"Only TP=$tpcounter% - $gtp% globally" "$1/finalsumresults.txt"

sed -i '1i\'"$fcnt	$fpcnt" "$1/finalsumresults.txt"


chmod go+rwx "$1/finalsumresults.txt"
chmod go+rwx "$1/finalsumtpresults.match.txt"
chmod go+rwx "$1/finalsumfpresults.match.txt"
chmod go+rwx "$1/finalsumtfpresults.match.txt"
chmod go+rwx "$1/finalsumresults.match.txt"