#!/bin/bash

BASEDIR=/home/mjwu/EteRNABot/benchmark/code
BIN=$BASEDIR/bin
TEMP=$BASEDIR/temp
RESFILE=$BASEDIR/results.txt
TIMEFILE=$BASEDIR/timing.txt

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

maxattempts=5

# RUN VIENNA
start=`date +%s`
i=0
solved=false
while [ $i -lt $maxattempts ]; do
    return=$(echo $secstruct | $BIN/vienna)
    if [[ $return == *"d="* ]]; then
        ((i++))
        continue
    else
        result+=$(echo $return | awk '{print $1}')"\t"
        solved=true
        break
    fi
done
if [ ! $solved ]; then
    result+="0\t"
fi
end=`date +%s`
time+=$((end-start))"\t"

# RUN INFORNA
start=`date +%s`
i=0
solved=false
while [ $i -lt $maxattempts ]; do
    return=$($BIN/inforna "$secstruct" 2> /dev/null | grep MFE)
    if [[ -z $return ]]; then
        ((i++))
        continue
#result+="0\t"
    elif [[ $return == *"NO_MFE"* ]]; then 
        ((i++))
        continue
#result+="0\t"
    else
    #    checkSolution $(echo $return | awk '{print $2}')
        result+=$(echo $return | awk '{print $2}')"\t"
        solved=true
        break
    fi
done
if [ ! $solved ]; then
    result+="0\t"
fi
end=`date +%s`
time+=$((end-start))"\t"

# RUN NUPACK
export NUPACKHOME=/home/mjwu/EteRNABot/benchmark/code/bin
start=`date +%s`
i=0
solved=false
while [ $i -lt $maxattempts ]; do
    echo $secstruct > $TEMP/${1}.fold
    #fstop=$(echo "scale=2;0.25/${#secstruct}" | bc)
    $BIN/nupack -material $BIN/rna1999 $TEMP/${1} > /dev/null
    violations=$(grep "Pattern" $TEMP/${1}.summary | awk '{print $4}')
    defect=$(grep "MFE defect" $TEMP/${1}.summary | awk '{print $4}')
    if [[ $violations == "0"  && $defect == "0.0"* ]]; then
        result+=$(tail -1 $TEMP/${1}.summary)"\t"
        solved=true
        break
    else
        ((i++))
        continue
    fi
done
if [ ! $solved ]; then
    result+="0\t"
fi
#rm $TEMP/${1}.summary
end=`date +%s`
time+=$((end-start))"\t"

# RUN RNA-SSD
start=`date +%s`
i=0
solved=false
while [ $i -lt $maxattempts ]; do
    return=$($BIN/RNAdesigner -f $TEMP/${1}.fold)
    dist=$(echo "$return" | grep "Distance" | awk '{print $2}')
    if [[ $dist == "0" ]]; then
        result+=$(echo "$return" | grep "Designed seq" | awk '{print $3}')"\t"
        solved=true
        break
    #    checkSolution $(echo "$return" | grep "Designed seq" | awk '{print $3}')
    else
        ((i++))
        continue
#result+="0\t"
    fi
done
if [ ! $solved ]; then
    result+="0\t"
fi
end=`date +%s`
time+=$((end-start))"\t"

# RUN MODENA
start=`date +%s`
i=0
solved=false
mkdir ${1}
cd ${1}
while [ $i -lt $maxattempts ]; do
    $BIN/modena -f $TEMP/${1}.fold > $TEMP/${1}.modena
    if grep -Fq "ERROR" $TEMP/${1}.modena; then
        break
    fi
    return=$(python $BASEDIR/parse_modena_out.py $TEMP/${1}.modena)
    if [[ -z $return ]]; then
        ((i++))
        continue
#result+="0\t"
    else
#checkSolution $return
        result+=$return"\t"
        solved=true
        break
    fi
done
if ! $solved; then
    result+="0\t"
fi
rm -rf ${1}
#rm $TEMP/${1}.fold
#rm $TEMP/${1}.modena
end=`date +%s`
time+=$((end-start))"\t"

# RUN DSS-OPT
start=`date +%s`
i=0
solved=false
while [ $i -lt $maxattempts ]; do
    $BIN/dss-opt $secstruct > $TEMP/${1}.dss
    return=$(python $BASEDIR/parse_dss_out.py $TEMP/${1}.dss)
    if [[ -z $return ]]; then
        ((i++))
        continue
#result+="0\t"
    else
    #    checkSolution $return
        result+="$return\t"
        solved=true
        break
    fi
done
if [ ! $solved ]; then
    result+="0\t"
fi
#rm $TEMP/${1}.dss
end=`date +%s`
time+=$((end-start))

echo -e $result >> $RESFILE
echo -e $time >> $TIMEFILE
echo "done"
