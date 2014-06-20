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

mkdir -p "$2/BOOST"
chmod go+rxw "$2/BOOST"
mkdir -p "$2/BEAM3"
chmod go+rxw "$2/BEAM3"
mkdir -p "$2/TEAM"
chmod go+rxw "$2/TEAM"
mkdir -p "$2/SNPR"
chmod go+rxw "$2/SNPR"
mkdir -p "$2/SNPH"
chmod go+rxw "$2/SNPH"

mkdir -p "$3/BOOST"
chmod go+rxw "$3/BOOST"
mkdir -p "$3/BEAM3"
chmod go+rxw "$3/BEAM3"
mkdir -p "$3/TEAM"
chmod go+rxw "$3/TEAM"
mkdir -p "$3/SNPR"
chmod go+rxw "$3/SNPR"
mkdir -p "$3/SNPH"
chmod go+rxw "$3/SNPH"

mkdir -p "$3/Ensemble"
chmod go+rxw "$3/Ensemble"

var1=$1
export var2=$2
export var3=$3
find $var1 -type f -name "*.mdr" | sort -n | 
xargs -i --max-procs=1 bash -c 'bash Ensembleformat.sh {} $var2 $var3;'