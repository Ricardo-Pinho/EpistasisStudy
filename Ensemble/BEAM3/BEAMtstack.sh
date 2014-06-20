#!/bin/bash



tcnt=0
cnt=0
cpu=0
tcpu=0
mem=0
tmem=0

> "$1/finaltresults.txt"

echo $tcnt > "$1/finaltresults.match.txt"
echo $cpu > "$1/finalcpuresults.match.txt"
echo $mem > "$1/finalmemresults.match.txt"

find $1 -type f -name "time*" |
while read file;
do

		echo ""	
		echo "***************************************************"
		echo "$file"

		while read line
		do

		if [[ $line == *"User time"* ]]; then

			cnt=${line:21:${#line}};
			tcnt=$(echo $tcnt + $cnt | bc);
			sed -i '1i\'"$file" "$1/finaltresults.match.txt"
			sed -i '1i\'"$tcnt" "$1/finaltresults.match.txt"

		elif [[ $line == *"Percent of CPU"* ]]; then
			
			cpu=${line:29:-1};
			tcpu=$(echo $tcpu + $cpu | bc);
			sed -i '1i\'"$file" "$1/finalcpuresults.match.txt"
			sed -i '1i\'"$tcpu" "$1/finalcpuresults.match.txt"
		elif [[ $line == *"Maximum resident set size"* ]]; then
			
			mem=${line:35:${#line}};
			tmem=$(echo $tmem + $mem | bc);
			sed -i '1i\'"$file" "$1/finalmemresults.match.txt"
			sed -i '1i\'"$tmem" "$1/finalmemresults.match.txt"
		fi	

		done < $file

		echo "Time: $cnt | CPU: $cpu% | Memory Usage: $mem (kbytes)"
		echo "***************************************************"
		echo ""

done > "$1/finaltresults.txt"


ttotal=$(sed -n 1p "$1/finaltresults.match.txt")
cputotal=$(sed -n 1p "$1/finalcpuresults.match.txt")
memtotal=$(sed -n 1p "$1/finalmemresults.match.txt")

ttavg=$(echo $ttotal*0.01 | bc)
tcpuavg=$(echo $cputotal*0.01 | bc)
tmemavg=$(echo $memtotal*0.01 | bc)

sed -i '1i\'"Average Time=$ttavg (seconds)" "$1/finaltresults.txt"
sed -i '1i\'"Average CPU Usage =$tcpuavg%" "$1/finaltresults.txt"
sed -i '1i\'"Average Memory Usage =$tmemavg (kbytes)" "$1/finaltresults.txt"

sed -i '1i\'"$ttavg $tcpuavg $tmemavg" "$1/finaltresults.txt"

chmod go+rwx "$1/finaltresults.txt"
chmod go+rwx "$1/finalcpuresults.match.txt"
chmod go+rwx "$1/finaltresults.match.txt"
chmod go+rwx "$1/finalmemresults.match.txt"