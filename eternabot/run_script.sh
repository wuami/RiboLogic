#!/bin/bash

#PBS -M mjwu@stanford.edu
#PBS -m ea

export PYTHONPATH=$PYTHONPATH:/home/mjwu/lib/python2.6/site-packages/
export PATH=$PATH:/home/mjwu/EteRNABot/eternabot/resources/vienna/linux

cd EteRNABot/eternabot
python switch_design.py and_gate
