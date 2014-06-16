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
export folderNo=70
export folderNo2=70
find $var1 -type f -name "*.mdr" | sort -n | xargs -i --max-procs=12 bash -c 'bash MBMDRformat.sh {} $format $var2 $folderNo $folderNo2; bash MBMDRconfig.sh {} $result $format $var2 $folderNo $folderNo2;'