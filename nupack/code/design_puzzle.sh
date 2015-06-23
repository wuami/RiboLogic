#!/bin/bash

id=${1}
OUTDIR=/home/mjwu/EteRNABot/nupack/data

if grep -Fxq $id $OUTDIR/done.txt
then
    echo "puzzle $id already done"
    exit 1
fi


echo "getting puzzle $id"
secstruct=$(curl -s "http://eternagame.org/get/?type=puzzle&nid=$id" | grep -E -o "secstruct.:.[\(\.\) ]+" | cut -c 13-)
secstruct=${secstruct//[[:blank:]]/}
locks=$(curl -s "http://eternagame.org/get/?type=puzzle&nid=$id" | grep -E -o "locks.:.[xo]+" | cut -c 9-)
if [ -z $locks ]; then
    n=${#secstruct}
    locks=$(printf "%0.so" $(seq 1 $n))
fi
printf "$secstruct\n$locks" > temp.fold


echo "designing rna"

$NUPACKHOME/bin/design -material rna temp;

python get_sequence_info.py temp.summary $id

#return=$(curl -s -X POST -c cookiefile -d "type=login&name=thenupackbot&pass=nilespierce&workbranch&main"  http://eternagame.org/login/)
#if [[ $return != *'"success":true'* ]]
#then
#    grep -P -o '"error":"(.*?)"' <<< $return
#fi
#return=$(curl -s -X POST -b cookiefile -d @temp.summary_post http://eternagame.org/post/)
#if [[ $return != *'"success":true'* ]]
#then
#    grep -P -o '"error":"(.*?)"' <<< $return
#fi

echo $id >> $OUTDIR/done.txt
