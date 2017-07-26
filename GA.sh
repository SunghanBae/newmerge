#!/bin/bash
#brips=1193
eurica=1408
aida=1208

for((i=0;i<8;i++))
do

name=GA_G${eurica}_A${aida}_${i}.root

if [ -e $name.root ]; then
        echo "$name.root exist"
else

time ./Beta_gamma $eurica $aida $i
fi

done

