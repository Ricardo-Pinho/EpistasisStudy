#!/bin/bash


if [[ $2 == *"/" ]]
then
   echo "error. Remove the '/' at the end."
   break
fi

if [[ $1 == *"/" ]]
then
   echo "error. Remove the '/' at the end."
   break
fi

if [[ $3 == *"/" ]]
then
   echo "error. Remove the '/' at the end."
   break
fi


var1=$1
export format=$2
export result=$3
export var2=$1
find $var1 -type f -name "*.mdr" | sort -n | 
xargs -i --max-procs=1 bash -c 'bash SnCformat.sh {} $format;
 bash SnCconfig.sh {} $result $format;'