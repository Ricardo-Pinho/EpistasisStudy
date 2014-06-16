#!/bin/bash

find $1 -type f -name "results*.txt.gz" |
while read file;
do

echo "$file"


gzip -dc "$file" > "${file%.gz}"

done 
