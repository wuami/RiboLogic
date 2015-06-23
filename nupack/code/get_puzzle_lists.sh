#!/bin/sh

#PBS -M mjwu@stanford.edu
#PBS -m ea
#PBS -l walltime=3:00:00:00

DATA_DIR=/home/mjwu/EteRNABot/nupack/data

curl "http://eternagame.org/get/?type=puzzles&sort=solved&puzzle_type=Challenge&skip=2&size=240&notcleared=true&uid=216510" -o $DATA_DIR/list.txt
grep -E -o "id.:.[0-9]+" $DATA_DIR/list.txt | cut -c 6- > $DATA_DIR/list.puzzles
curl "http://eternagame.org/get/?type=puzzles&sort=date&puzzle_type=PlayerPuzzle&skip=0&size=8000&notcleared=true&uid=216510" -o $DATA_DIR/list.txt
grep -E -o "id.:.[0-9]+" $DATA_DIR/list.txt | cut -c 6- > $DATA_DIR/list.player
curl "http://eternagame.org/get/?type=puzzles&sort=solved&puzzle_type=PlayerPuzzle&skip=0&size=8000&notcleared=true&uid=216510" -o $DATA_DIR/list.txt
grep -E -o "id.:.[0-9]+" $DATA_DIR/list.txt | cut -c 6- > $DATA_DIR/list.easy

xargs < $DATA_DIR/list.puzzles -I % ./design_puzzle.sh %
xargs < $DATA_DIR/list.player -I % ./design_puzzle.sh %
xargs < $DATA_DIR/list.easy -I % ./design_puzzle.sh %
