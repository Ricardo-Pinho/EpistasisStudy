#!/bin/bash


index=1;
chgfile="test";
file=$1

newfile="$2/$(basename $(dirname ${file}))/$(basename ${file})"

if [[ ! -f "${newfile}.txt" ]]
then

if [[ ! $file == *$chgfile* ]];
then
	mkdir -p "$2/$(basename $(dirname ${file}))"
	chgfile="$(dirname ${file})";
	chmod go+rxw "$2/$(basename $(dirname ${file}))"
fi

index=$((index+1))
echo "$index";

tr -s '[:blank:]' ' ' < ${file} > ${newfile}.temp.txt;

sed -i '1i\'"Pos 1 2 3 4 5 6 7 8 9 10 11 12 13 14 15 16 17 18 19 20 21 22 23 24 25 26 27 28 29 30 31 32 33 34 35 36 37 38 39 40 41 42 43 44 45 46 47 48 49 50 51 52 53 54 55 56 57 58 59 60 61 62 63 64 65 66 67 68 69 70 71 72 73 74 75 76 77 78 79 80 81 82 83 84 85 86 87 88 89 90 91 92 93 94 95 96 97 98 99 100 101 102 103 104 105 106 107 108 109 110 111 112 113 114 115 116 117 118 119 120 121 122 123 124 125 126 127 128 129 130 131 132 133 134 135 136 137 138 139 140 141 142 143 144 145 146 147 148 149 150 151 152 153 154 155 156 157 158 159 160 161 162 163 164 165 166 167 168 169 170 171 172 173 174 175 176 177 178 179 180 181 182 183 184 185 186 187 188 189 190 191 192 193 194 195 196 197 198 199 200 201 202 203 204 205 206 207 208 209 210 211 212 213 214 215 216 217 218 219 220 221 222 223 224 225 226 227 228 229 230 231 232 233 234 235 236 237 238 239 240 241 242 243 244 245 246 247 248 249 250 251 252 253 254 255 256 257 258 259 260 261 262 263 264 265 266 267 268 269 270 271 272 273 274 275 276 277 278 279 280 281 282 283 284 285 286 287 288 289 290 291 292 293 294 295 296 297 298 299 300" ${newfile}.temp.txt

sed -i '1i\'"Chr chr1 chr1 chr1 chr1 chr1 chr1 chr1 chr1 chr1 chr1 chr1 chr1 chr1 chr1 chr1 chr1 chr1 chr1 chr1 chr1 chr1 chr1 chr1 chr1 chr1 chr1 chr1 chr1 chr1 chr1 chr1 chr1 chr1 chr1 chr1 chr1 chr1 chr1 chr1 chr1 chr1 chr1 chr1 chr1 chr1 chr1 chr1 chr1 chr1 chr1 chr1 chr1 chr1 chr1 chr1 chr1 chr1 chr1 chr1 chr1 chr1 chr1 chr1 chr1 chr1 chr1 chr1 chr1 chr1 chr1 chr1 chr1 chr1 chr1 chr1 chr1 chr1 chr1 chr1 chr1 chr1 chr1 chr1 chr1 chr1 chr1 chr1 chr1 chr1 chr1 chr1 chr1 chr1 chr1 chr1 chr1 chr1 chr1 chr1 chr1 chr1 chr1 chr1 chr1 chr1 chr1 chr1 chr1 chr1 chr1 chr1 chr1 chr1 chr1 chr1 chr1 chr1 chr1 chr1 chr1 chr1 chr1 chr1 chr1 chr1 chr1 chr1 chr1 chr1 chr1 chr1 chr1 chr1 chr1 chr1 chr1 chr1 chr1 chr1 chr1 chr1 chr1 chr1 chr1 chr1 chr1 chr1 chr1 chr1 chr1 chr1 chr1 chr1 chr1 chr1 chr1 chr1 chr1 chr1 chr1 chr1 chr1 chr1 chr1 chr1 chr1 chr1 chr1 chr1 chr1 chr1 chr1 chr1 chr1 chr1 chr1 chr1 chr1 chr1 chr1 chr1 chr1 chr1 chr1 chr1 chr1 chr1 chr1 chr1 chr1 chr1 chr1 chr1 chr1 chr1 chr1 chr1 chr1 chr1 chr1 chr2 chr2 chr2 chr2 chr2 chr2 chr2 chr2 chr2 chr2 chr2 chr2 chr2 chr2 chr2 chr2 chr2 chr2 chr2 chr2 chr2 chr2 chr2 chr2 chr2 chr2 chr2 chr2 chr2 chr2 chr2 chr2 chr2 chr2 chr2 chr2 chr2 chr2 chr2 chr2 chr2 chr2 chr2 chr2 chr2 chr2 chr2 chr2 chr2 chr2 chr2 chr2 chr2 chr2 chr2 chr2 chr2 chr2 chr2 chr2 chr2 chr2 chr2 chr2 chr2 chr2 chr2 chr2 chr2 chr2 chr2 chr2 chr2 chr2 chr2 chr2 chr2 chr2 chr2 chr2 chr2 chr2 chr2 chr2 chr2 chr2 chr2 chr2 chr2 chr2 chr2 chr2 chr2 chr2 chr2 chr2 chr2 chr2 chr2 chr2" ${newfile}.temp.txt

sed -i '1i\'"ID rs1 rs2 rs3 rs4 rs5 rs6 rs7 rs8 rs9 rs10 rs11 rs12 rs13 rs14 rs15 rs16 rs17 rs18 rs19 rs20 rs21 rs22 rs23 rs24 rs25 rs26 rs27 rs28 rs29 rs30 rs31 rs32 rs33 rs34 rs35 rs36 rs37 rs38 rs39 rs40 rs41 rs42 rs43 rs44 rs45 rs46 rs47 rs48 rs49 rs50 rs51 rs52 rs53 rs54 rs55 rs56 rs57 rs58 rs59 rs60 rs61 rs62 rs63 rs64 rs65 rs66 rs67 rs68 rs69 rs70 rs71 rs72 rs73 rs74 rs75 rs76 rs77 rs78 rs79 rs80 rs81 rs82 rs83 rs84 rs85 rs86 rs87 rs88 rs89 rs90 rs91 rs92 rs93 rs94 rs95 rs96 rs97 rs98 rs99 rs100 rs101 rs102 rs103 rs104 rs105 rs106 rs107 rs108 rs109 rs110 rs111 rs112 rs113 rs114 rs115 rs116 rs117 rs118 rs119 rs120 rs121 rs122 rs123 rs124 rs125 rs126 rs127 rs128 rs129 rs130 rs131 rs132 rs133 rs134 rs135 rs136 rs137 rs138 rs139 rs140 rs141 rs142 rs143 rs144 rs145 rs146 rs147 rs148 rs149 rs150 rs151 rs152 rs153 rs154 rs155 rs156 rs157 rs158 rs159 rs160 rs161 rs162 rs163 rs164 rs165 rs166 rs167 rs168 rs169 rs170 rs171 rs172 rs173 rs174 rs175 rs176 rs177 rs178 rs179 rs180 rs181 rs182 rs183 rs184 rs185 rs186 rs187 rs188 rs189 rs190 rs191 rs192 rs193 rs194 rs195 rs196 rs197 rs198 rs199 rs200 rs201 rs202 rs203 rs204 rs205 rs206 rs207 rs208 rs209 rs210 rs211 rs212 rs213 rs214 rs215 rs216 rs217 rs218 rs219 rs220 rs221 rs222 rs223 rs224 rs225 rs226 rs227 rs228 rs229 rs230 rs231 rs232 rs233 rs234 rs235 rs236 rs237 rs238 rs239 rs240 rs241 rs242 rs243 rs244 rs245 rs246 rs247 rs248 rs249 rs250 rs251 rs252 rs253 rs254 rs255 rs256 rs257 rs258 rs259 rs260 rs261 rs262 rs263 rs264 rs265 rs266 rs267 rs268 rs269 rs270 rs271 rs272 rs273 rs274 rs275 rs276 rs277 rs278 rs279 rs280 rs281 rs282 rs283 rs284 rs285 rs286 rs287 rs288 rs289 rs290 rs291 rs292 rs293 rs294 rs295 rs296 rs297 rs298 rs299 rs300" ${newfile}.temp.txt

./transpose.awk "${newfile}.temp.txt" > "${newfile}.txt"

rm ${newfile}.temp.txt;

chmod go+rwx ${newfile}.txt

fi