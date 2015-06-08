import design_sequence
from design_utils import bp_distance_with_unpaired as bp

f = open("resources/puzzles/xor_gate_3parity1.json", 'r')
p = design_sequence.read_puzzle_json(f.read())
s1 = ".....................................................(((((((((((((((.....................................................................)))))))...))))...))))"
print bp(s1, p.targets[0]['secstruct'], p.targets[0]['constrained'], p.targets[0]['threshold'])
