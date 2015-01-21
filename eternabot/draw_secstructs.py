import argparse
import settings
import os
import varna
from switch_design import read_puzzle_json
import inv_utils

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
colors = 1.0/(1+n)
colormaps = []
for target in puzzle.targets:
     colormap = "1.0;"*puzzle.n
     for i, key in enumerate(sorted(inputs)):
         if i != 0:
            colormap += "0.0;"*puzzle.linker_length
         if key in target['inputs']:
             colormap += (str((i+1)*colors)+';')*len(inputs[key])
         else:
             colormap += "0.0;"*len(inputs[key])
     print colormap
     colormaps.append(colormap)

# draw image for each sequence
n_sequences = 0
with open(os.path.join(settings.PUZZLE_DIR, "%s.out" % args.puzzleid), 'r') as f:
    for line in f:
        if not line.startswith("#"):
            seq = line.split()[0]
            for j, target in enumerate(puzzle.targets):
                foldseq = puzzle.get_fold_sequence(seq, target)
                secstruct = inv_utils.fold(foldseq)[0]
                foldseq = foldseq.replace("&", "\&")
                filename = "%s/images/%s_%s-%s.png" % (settings.PUZZLE_DIR, args.puzzleid, n_sequences, j)
                v.new_image_by_str(filename, secstruct, foldseq, colormap_str=colormaps[j])
            n_sequences += 1

# create html to display images
htmlfile = os.path.join(settings.PUZZLE_DIR, "images", "%s.html" % args.puzzleid)
with open(htmlfile, 'w') as f: 
    f.write("<html>\n<body>\n")
    f.write("<table>\n")
    for i in range(n_sequences):
        f.write("<tr>\n")
        for j in range(puzzle.n_targets):
            f.write("\t<td><img src=\"%s_%s-%s.png\" style=\"max-width:500px;max-height:500px\"></td>\n" % (args.puzzleid, i, j))
        f.write("</tr>\n")
    f.write("</body>\n</html>")
