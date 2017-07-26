#!/bin/bash
brips=1193
#eurica=1408
aida=1208

for((i=0;i<8;i++))
do

name=BA_B${brips}_A${aida}_${i}.root

if [ -e $name.root ]; then
        echo "$name.root exist"
else

time ./Beam_Ion $brips $aida $i
fi

done

