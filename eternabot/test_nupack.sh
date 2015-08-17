#PBS -M mjwu@stanford.edu
#PBS -m ea
#PBS -l walltime=2:00:00:00
#PBS -l nodes=1:ppn=1

source ~/.bashrc
python2.7 ~/EteRNABot/eternabot/design_sequence.py -m nupack and_gate_new
