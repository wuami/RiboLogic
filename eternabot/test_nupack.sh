#!/bin/bash

#PBS -M mjwu@stanford.edu
#PBS -m ea
#PBS -l walltime=2:00:00:00
#PBS -l nodes=1:ppn=4

source ~/.bashrc

puzzle="and_gate_tb"
niter=10000
echo $puzzle", $niter iter"
python2.7 /home/mjwu/EteRNABot/eternabot/design_sequence.py -m nupack --prints -i $niter $puzzle
