#!/bin/bash

awk '/================/{p=0}p;/======Significant SNP-pairs=====/{p=1}' "/home/ricardo/Desktop/TEAMResults/0.5,2000,I,2.0,0.02/results48.txt" > "/home/ricardo/Desktop/TEAMResults/0.5,2000,I,2.0,0.02/temp.txt"

while read line
		do
			echo  "${line#*:} ${line%:*}"
		done < "/home/ricardo/Desktop/TEAMResults/0.5,2000,I,2.0,0.02/temp.txt" > "/home/ricardo/Desktop/TEAMResults/0.5,2000,I,2.0,0.02/temp2.txt"

		sort -r -g "/home/ricardo/Desktop/TEAMResults/0.5,2000,I,2.0,0.02/temp2.txt" > "/home/ricardo/Desktop/TEAMResults/0.5,2000,I,2.0,0.02/temp.txt"
