#!/bin/bash

BASEDIR=/home/mjwu/EteRNABot/benchmark
CODEDIR=$BASEDIR/code
DATADIR=$BASEDIR/data
BIN=$CODEDIR/bin
TEMP=$DATADIR/temp
RESFILE=$DATADIR/results_dss.txt
TIMEFILE=$DATADIR/timing_dss.txt

source ~/.bashrc

#if grep -Fq ${1} $RESFILE; then
#    echo "${1} already run"
#    exit 1
#fi
echo -n "running ${1}..."

secstruct=$(grep -F ${1} $DATADIR/EteRNA100.txt | awk '{print $3}')
result=${1}"\t"
time=${1}"\t"${#secstruct}"\t"

function checkSolution {
    sequence=$1
    fold=$(echo "$sequence" | $BIN/RNAfold | tail -1 | awk '{print $1}')
    if [[ $fold == $secstruct ]]; then
        result+="$sequence\t"
        solved=true
    fi
}

maxattempts=5

# RUN DSS-OPT
start=`date +%s`
i=0
solved=false
while [ $i -lt $maxattempts ]; do
    echo $i
    $BIN/dss-opt $secstruct > $TEMP/${1}.dss
    return=$(python $CODEDIR/parse_dss_out.py $TEMP/${1}.dss)
    if [[ -z $return ]]; then
        ((i++))
        continue
#result+="0\t"
    else
        checkSolution $return
        if $solved; then
            break
        else
            ((i++))
            continue
        fi
    fi
done
if ! $solved; then
    result+="0\t"
fi
#rm $TEMP/${1}.dss
end=`date +%s`
time+=$((end-start))

echo -e $result >> $RESFILE
echo -e $time >> $TIMEFILE
echo "done"
