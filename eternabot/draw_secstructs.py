import argparse
import settings
import os
import varna
from design_sequence import read_puzzle_json
import inv_utils

def get_colormaps(targets, inputs, sequence_length, linker_length, n):
    if n == 2:
        colors = [0.5, 0.67]
    else:
        colors = [(i+1)/(1+n) for i in range(n)]
    colors = [str(x)+";" for x in colors]
    colormaps = []
    for target in targets:
         colormap = "1.0;"*sequence_length
         for i, key in enumerate(sorted(inputs)):
             if i != 0:
                colormap += "0.0;"*linker_length
             if key in target['inputs']:
                 colormap += colors[i]*len(inputs[key])
             else:
                 colormap += "0.0;"*len(inputs[key])
         colormaps.append(colormap)
    return colormaps

def draw_secstruct_state(v, target, foldseq, colormap, filename):
    secstruct = inv_utils.fold(foldseq)[0]
    foldseq = foldseq.replace("&", "\&")
    highlight = ""
    if secstruct[39:58] == "(((((.((....)))))))":
        highlight = "40-58:fill=#666666,outline=#FFFFFF,radius=20.0"
    v.new_image_by_str(filename, secstruct, foldseq, highlight_region=highlight, colormap_str=colormap)

def write_html(puzzleid, n_sequences, n_targets, order=None, breaks=[]):
    htmlfile = os.path.join(settings.PUZZLE_DIR, "images", "%s.html" % puzzleid)
    if order == None:
        order = range(n_sequences)
    
    with open(htmlfile, 'w') as f: 
        f.write("<html>\n<body>\n")
        f.write("<table rules=\"rows\">\n")
        f.write("<font size=\"20\">")
        f.write("<tr><td align=\"center\">no oligos</td><td align=\"center\">oligo 1</td><td align=\"center\">oligo 2</td><td align=\"center\">oligo 1 + oligo 2</td></tr>")
        f.write("</font>")
        for i,n in enumerate(order):
            f.write("<tr>\n")
            for j in range(n_targets):
                f.write("\t<td><img src=\"%s_%s-%s.png\" style=\"max-width:500px;max-height:500px\"></td>\n" % (puzzleid, n, j))
            f.write("</tr>\n")
            if i in breaks:
                f.write("</table><br>\n")
                f.write("<hr style=\"color:black;height:10px\"/>")
                f.write("Cluster %s\n" % breaks.index(i))
                f.write("<table rules=\"rows\">\n")
        f.write("</table>\n</body>\n</html>")


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
    colormaps = get_colormaps(puzzle.targets, inputs, puzzle.n, puzzle.linker_length, n)
    
    # draw image for each sequence
    n_sequences = 0
    with open(os.path.join(settings.PUZZLE_DIR, "%s.out" % args.puzzleid), 'r') as f:
        for line in f:
            if not line.startswith("#"):
                seq = line.split()[0]
                for j, target in enumerate(puzzle.targets):
                    filename = "%s/images/%s_%s-%s.png" % (settings.PUZZLE_DIR, args.puzzleid, n_sequences, j)
                    draw_secstruct_state(v, target, puzzle.get_fold_sequence(seq, target), colormaps[j], filename)
                n_sequences += 1
    
    # create html to display images
    write_html(args.puzzleid, n_sequences, n_targets)

if __name__ == "__main__":
    main()
