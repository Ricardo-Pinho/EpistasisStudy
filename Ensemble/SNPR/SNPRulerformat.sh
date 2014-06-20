#!/bin/bash


index=1;
chgfile="test";
file=$1

newfile="$2/$(basename $(dirname ${file}))/$(basename ${file})"

temp=$(echo "$file" | sed -e 's/.*1750.1.//g')
ind=$(echo "$temp" | sed -e 's/-3GE-.*//g')
resultfolder="$3/$(basename $(dirname ${file}))"

if [[ ! -f "${newfile}.txt" ]] && [[ ! -f "$resultfolder/results${ind}.txt" ]]
then

echo $file
echo $chgfile

if [[ ! $file == *$chgfile* ]];
then
	mkdir -p "$2/$(basename $(dirname ${file}))"
	chgfile="$(dirname ${file})";
	chmod go+rxw "$2/$(basename $(dirname ${file}))"
fi

while read a;
do
fchar=${a:0:1};
#line = ${a:0:${#a}}
linex=${a:4:${#a}};
echo "$linex $fchar";
done < $file > ${newfile}.temp.txt

index=$((index+1))
echo "$index";

tr -s '[:blank:]' ' ' < ${newfile}.temp.txt > ${newfile}.txt;
sed -i '1i\'"X1 X2 X3 X4 X5 X6 X7 X8 X9 X10 X11 X12 X13 X14 X15 X16 X17 X18 X19 X20 X21 X22 X23 X24 X25 X26 X27 X28 X29 X30 X31 X32 X33 X34 X35 X36 X37 X38 X39 X40 X41 X42 X43 X44 X45 X46 X47 X48 X49 X50 X51 X52 X53 X54 X55 X56 X57 X58 X59 X60 X61 X62 X63 X64 X65 X66 X67 X68 X69 X70 X71 X72 X73 X74 X75 X76 X77 X78 X79 X80 X81 X82 X83 X84 X85 X86 X87 X88 X89 X90 X91 X92 X93 X94 X95 X96 X97 X98 X99 X100 X101 X102 X103 X104 X105 X106 X107 X108 X109 X110 X111 X112 X113 X114 X115 X116 X117 X118 X119 X120 X121 X122 X123 X124 X125 X126 X127 X128 X129 X130 X131 X132 X133 X134 X135 X136 X137 X138 X139 X140 X141 X142 X143 X144 X145 X146 X147 X148 X149 X150 X151 X152 X153 X154 X155 X156 X157 X158 X159 X160 X161 X162 X163 X164 X165 X166 X167 X168 X169 X170 X171 X172 X173 X174 X175 X176 X177 X178 X179 X180 X181 X182 X183 X184 X185 X186 X187 X188 X189 X190 X191 X192 X193 X194 X195 X196 X197 X198 X199 X200 X201 X202 X203 X204 X205 X206 X207 X208 X209 X210 X211 X212 X213 X214 X215 X216 X217 X218 X219 X220 X221 X222 X223 X224 X225 X226 X227 X228 X229 X230 X231 X232 X233 X234 X235 X236 X237 X238 X239 X240 X241 X242 X243 X244 X245 X246 X247 X248 X249 X250 X251 X252 X253 X254 X255 X256 X257 X258 X259 X260 X261 X262 X263 X264 X265 X266 X267 X268 X269 X270 X271 X272 X273 X274 X275 X276 X277 X278 X279 X280 X281 X282 X283 X284 X285 X286 X287 X288 X289 X290 X291 X292 X293 X294 X295 X296 X297 X298 X299 X300 Label" ${newfile}.txt
rm ${newfile}.temp.txt;
chmod go+rwx ${newfile}.txt
fi