#!/bin/bash

BASEDIR=/home/mjwu/EteRNABot/benchmark/code
BIN=$BASEDIR/bin
TEMP=$BASEDIR/temp
RESFILE=$BASEDIR/results_modena.txt
TIMEFILE=$BASEDIR/timing_modena.txt

source ~/.bashrc

if grep -Fq ${1} $RESFILE; then
    echo "${1} already run"
    exit 1
fi
echo -n "running ${1}..."

secstruct=$(grep -F ${1} $BASEDIR/puzzles.txt | awk '{print $3}')
result=${1}"\t"
time=${1}"\t"${#secstruct}"\t"

function checkSolution {
    sequence=$1
    fold=$(echo "$sequence" | $BIN/RNAfold | tail -1 | awk '{print $1}')
    if [[ $fold == $secstruct ]]; then
        result+="$sequence\t"
    else
        result+="0\t"
    fi
}

# RUN MODENA
start=`date +%s`
$BIN/modena -f $TEMP/${1}.fold > $TEMP/${1}.modena
return=$(python $BASEDIR/parse_modena_out.py $TEMP/${1}.modena)
if [[ -z $return ]]; then
    result+="0\t"
else
    checkSolution $return
#    result+=$return"\t"
fi
#rm $TEMP/${1}.fold
#rm $TEMP/${1}.modena
end=`date +%s`
time+=$((end-start))"\t"

echo -e $result >> $RESFILE
echo -e $time >> $TIMEFILE
echo "done"
