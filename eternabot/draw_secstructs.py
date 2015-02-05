import argparse
import varna
import settings, os
from design_sequence import read_puzzle_json
import draw_utils

def main():
    # parse arguments
    p = argparse.ArgumentParser()
    p.add_argument('puzzleid', help="name of puzzle filename", type=str)
    args = p.parse_args()
    
    # initiate varna RNA visualizer
    v = varna.Varna()
    
    # get puzzle object and generate colormaps for each objective
    json = open(os.path.join(settings.PUZZLE_DIR, "%s.json" % args.puzzleid)).read()
    puzzle = read_puzzle_json(json)
    inputs = puzzle.inputs
    n = len(inputs)
    colormaps = draw_utils.get_colormaps(puzzle.targets, inputs, puzzle.n, puzzle.linker_length, puzzle.design_linker, n)
    
    # draw image for each sequence
    n_sequences = 0
    with open(os.path.join(settings.PUZZLE_DIR, "%s.out" % args.puzzleid), 'r') as f:
        for line in f:
            if not line.startswith("#"):
                seq = line.split()[0]
                for j, target in enumerate(puzzle.targets):
                    filename = "%s/images/%s_%s-%s.png" % (settings.PUZZLE_DIR, args.puzzleid, n_sequences, j)
                    draw_utils.draw_secstruct_state(v, target, puzzle.get_fold_sequence(seq, target), colormaps[j], filename)
                n_sequences += 1
    
    # create html to display images
    draw_utils.write_html(args.puzzleid, n_sequences, n_targets)

if __name__ == "__main__":
    main()
