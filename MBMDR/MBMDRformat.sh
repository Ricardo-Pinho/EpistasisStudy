#!/bin/bash


index=1;
chgfile="test";
file=$1
folderPath=$3;
#folderNo=$4;
#folderNo2=$5;
newfile="$2/$(basename $(dirname ${file}))/$(basename ${file})"


if [[ ! $file == *$chgfile* ]];
then
	mkdir -p "$2/$(basename $(dirname ${file}))"
	chgfile="$(dirname ${file})";
	chmod go+rxw "$2/$(basename $(dirname ${file}))"
fi

#directory=$( ls -1 $folderPath | sort -n | head -n $folderNo2 | tail -n $folderNo )

if [[ ! -f "${newfile}.csv" ]] #&& [[ $directory == *$(basename $(dirname ${file}))*   ]]
then

echo $file

lindex=0;
while read a;
do
echo "\"$lindex\" $a";
lindex=$((lindex+1));
done < $file > ${newfile}.temp.txt


tr -s '[:blank:]' ',' < ${newfile}.temp.txt > ${newfile}.csv;

sed -i '1i\''"","Y","SNP1","SNP2","SNP3","SNP4","SNP5","SNP6","SNP7","SNP8","SNP9","SNP10","SNP11","SNP12","SNP13","SNP14","SNP15","SNP16","SNP17","SNP18","SNP19","SNP20","SNP21","SNP22","SNP23","SNP24","SNP25","SNP26","SNP27","SNP28","SNP29","SNP30","SNP31","SNP32","SNP33","SNP34","SNP35","SNP36","SNP37","SNP38","SNP39","SNP40","SNP41","SNP42","SNP43","SNP44","SNP45","SNP46","SNP47","SNP48","SNP49","SNP50","SNP51","SNP52","SNP53","SNP54","SNP55","SNP56","SNP57","SNP58","SNP59","SNP60","SNP61","SNP62","SNP63","SNP64","SNP65","SNP66","SNP67","SNP68","SNP69","SNP70","SNP71","SNP72","SNP73","SNP74","SNP75","SNP76","SNP77","SNP78","SNP79","SNP80","SNP81","SNP82","SNP83","SNP84","SNP85","SNP86","SNP87","SNP88","SNP89","SNP90","SNP91","SNP92","SNP93","SNP94","SNP95","SNP96","SNP97","SNP98","SNP99","SNP100","SNP101","SNP102","SNP103","SNP104","SNP105","SNP106","SNP107","SNP108","SNP109","SNP110","SNP111","SNP112","SNP113","SNP114","SNP115","SNP116","SNP117","SNP118","SNP119","SNP120","SNP121","SNP122","SNP123","SNP124","SNP125","SNP126","SNP127","SNP128","SNP129","SNP130","SNP131","SNP132","SNP133","SNP134","SNP135","SNP136","SNP137","SNP138","SNP139","SNP140","SNP141","SNP142","SNP143","SNP144","SNP145","SNP146","SNP147","SNP148","SNP149","SNP150","SNP151","SNP152","SNP153","SNP154","SNP155","SNP156","SNP157","SNP158","SNP159","SNP160","SNP161","SNP162","SNP163","SNP164","SNP165","SNP166","SNP167","SNP168","SNP169","SNP170","SNP171","SNP172","SNP173","SNP174","SNP175","SNP176","SNP177","SNP178","SNP179","SNP180","SNP181","SNP182","SNP183","SNP184","SNP185","SNP186","SNP187","SNP188","SNP189","SNP190","SNP191","SNP192","SNP193","SNP194","SNP195","SNP196","SNP197","SNP198","SNP199","SNP200","SNP201","SNP202","SNP203","SNP204","SNP205","SNP206","SNP207","SNP208","SNP209","SNP210","SNP211","SNP212","SNP213","SNP214","SNP215","SNP216","SNP217","SNP218","SNP219","SNP220","SNP221","SNP222","SNP223","SNP224","SNP225","SNP226","SNP227","SNP228","SNP229","SNP230","SNP231","SNP232","SNP233","SNP234","SNP235","SNP236","SNP237","SNP238","SNP239","SNP240","SNP241","SNP242","SNP243","SNP244","SNP245","SNP246","SNP247","SNP248","SNP249","SNP250","SNP251","SNP252","SNP253","SNP254","SNP255","SNP256","SNP257","SNP258","SNP259","SNP260","SNP261","SNP262","SNP263","SNP264","SNP265","SNP266","SNP267","SNP268","SNP269","SNP270","SNP271","SNP272","SNP273","SNP274","SNP275","SNP276","SNP277","SNP278","SNP279","SNP280","SNP281","SNP282","SNP283","SNP284","SNP285","SNP286","SNP287","SNP288","SNP289","SNP290","SNP291","SNP292","SNP293","SNP294","SNP295","SNP296","SNP297","SNP298","SNP299","SNP300"' ${newfile}.csv

rm ${newfile}.temp.txt;

chmod go+rwx ${newfile}.csv
fi