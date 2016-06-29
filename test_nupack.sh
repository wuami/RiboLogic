#!/bin/bash

#PBS -M mkap@stanford.edu
#PBS -m ea
#PBS -l walltime=2:00:00:00
#PBS -l nodes=1:ppn=4

source ~/.bashrc

niter=100000
conc="1"
echo ${puzzle}", $niter iter, start conc $conc, BPP SCORING"
echo ${options}
python2.7 /home/mkap/Das/design_sequence.py -m nupack --print_ -i $niter -c $conc ${options} ${puzzle}

