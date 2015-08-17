#!/bin/bash

#PBS -M mjwu@stanford.edu
#PBS -m ea
#PBS -l walltime=2:00:00:00
#PBS -l nodes=1:ppn=16

CODEDIR=/home/mjwu/EteRNABot/benchmark/code
DATADIR=/home/mjwu/EteRNABot/benchmark/data

#echo "run algs1 fixed solve"
echo "run 5 timeout"
#xargs < $CODEDIR/puzzle_names.txt -P 16 -I % $CODEDIR/run_algs2.sh "%"
xargs < $DATADIR/EteRNA100_names.txt -P 16 -I % $CODEDIR/run_algs2.sh "%"

#echo "run notnupack"
#xargs < $DATADIR/todo.txt -P 16 -I % $CODEDIR/run_notnupack.sh "%"
