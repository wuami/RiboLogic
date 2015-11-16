import fold_utils

def get_colormaps(targets, inputs, sequence_length, n):
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
        colormap += colors[i]*len(inputs[key])
    colormap += "1.0;"*(sequence_length)
    #colormaps.append(colormap)
    return colormap

def draw_secstruct_state(v, target, foldseq, colormap, filename):
    if "fold_constraint" in target:
        secstruct = fold_utils.fold(foldseq, target['fold_constraint'])[0]
    else:
        secstruct = fold_utils.fold(foldseq)[0]
    foldseq = foldseq.replace("&", "\&")
    v.new_image_by_str(filename, secstruct, foldseq, colormap_str=colormap)
