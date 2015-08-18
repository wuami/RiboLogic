#!/bin/bash

sequence=$(cat /dev/urandom | tr -dc 'AUGC' | fold -w 50 | head -n 1)

vienna=$(echo $sequence | ../../benchmark/code/bin/RNAfold185 | tail -1)
[[ "$vienna" =~ ([\(\)\.]+)[[:space:]]\([[:space:]]*([-\.0-9]+)\) ]] #\(([-\.0-9]+)\) ]] #\s\(([0-9.]+)\)" ]]
vienna_fold=${BASH_REMATCH[1]}
vienna_energy=${BASH_REMATCH[2]}

echo $sequence > temp.in
nupack=$(nupack3.0.5c-clean/bin/mfe -material rna1999 temp)
nupack_fold=$(sed '16q;d' temp.mfe)
nupack_energy=$(sed '15q;d' temp.mfe)
#rm temp.in temp.mfe

if [ "$vienna_fold" = "$nupack_fold" ]; then # && $vienna_energy -eq $nupack_energy ]]; then
    if [[ ! $(echo "$vienna_energy == $nupack_energy" | bc) ]]; then
        echo "energy different"
        echo "vienna: "$vienna_energy
        echo "nupack: "$nupack_energy
    fi
else
    echo "fold different"
    echo "vienna: $vienna_fold $vienna_energy" 
    echo "nupack: $nupack_fold $nupack_energy"
fi
