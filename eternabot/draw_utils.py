import argparse
import settings
import os
import varna
import fold_utils

def get_colormaps(targets, inputs, input_pos, sequence_length, n):
    if n == 2:
        colors = [0.5, 0.65]
    elif n == 3:
        colors = [0.5, 0.65, 0.8]
    else:
        colors = [float(i)/(2*n)+0.4 for i in range(n)]
    colors = [str(x)+";" for x in colors]
    colormaps = []
    #for target in targets:
    colormap = ""
    for i, key in enumerate(sorted(inputs)):
        #if key in target['inputs']:
        #    colormap += colors[i]*len(inputs[key])
        #else:
        #    colormap += "0.0;"*len(inputs[key])
        colormap += "1.0;"*(input_pos[i+1]-input_pos[i])
        colormap += colors[i]*len(inputs[key])
    colormap += "1.0;"*(sequence_length-input_pos[-1])
    #colormaps.append(colormap)
    return colormap

def draw_secstruct_state(v, target, foldseq, colormap, filename):
    if "fold_constraint" in target:
        secstruct = fold_utils.fold(foldseq, target['fold_constraint'])[0]
    else:
        secstruct = fold_utils.fold(foldseq)[0]
    print secstruct
    foldseq = foldseq.replace("&", "\&")
    highlight = ""
    idx = secstruct.find("(((((.((....)))))))")
    if idx != -1:
        lo, hi = idx-1, idx+20
        highlight = "%s-%s:fill=#666666,outline=#FFFFFF,radius=20.0" % (lo+1, hi)
    v.new_image_by_str(filename, secstruct, foldseq, highlight_region=highlight, colormap_str=colormap)

def write_html(puzzleid, n_sequences, n_targets, design_order=None, state_order=None, breaks=[]):
    htmlfile = os.path.join(settings.PUZZLE_DIR, "images", "%s.html" % puzzleid)
    if design_order == None:
        design_order = range(n_sequences)
    if state_order == None:
        state_order = range(n_targets)
    if n_targets == 4:
        state_names = ["no oligos", "oligo 1", "oligo 2", "oligo 1 + oligo 2"]
    else:
        state_names = range(n_targets)
    
    with open(htmlfile, 'w') as f: 
        f.write("<html>\n<body>\n")
        f.write("<table rules=\"rows\" style=\"font-size:40pt\">\n")
        f.write("<tr><td></td>")
        for i in state_order:
            f.write("<td align=\"center\">%s</td>" % state_names[i])
        f.write("</tr>")
        for i,n in enumerate(design_order):
            f.write("<tr><td>%s</td>\n" % i)
            for j,m in enumerate(state_order):
                f.write("\t<td><img src=\"%s_%s-%s.png\" style=\"max-width:500px;max-height:500px\"></td>\n" % (puzzleid, n, m))
            f.write("</tr>\n")
            if i in breaks:
                f.write("</table><br>\n")
                f.write("<hr style=\"color:black;height:10px\"/>")
                f.write("Cluster %s\n" % breaks.index(i))
                f.write("<table rules=\"rows\">\n")
        f.write("</table>\n</body>\n</html>")

