#PBS -M mjwu@stanford.edu
#PBS -m ea
#PBS -l walltime=2:00:00:00
#PBS -l nodes=1:ppn=16

source ~/.bashrc
python2.7 ~/EteRNABot/eternabot/design_sequence.py -m nupack xor_gate_new
