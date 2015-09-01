import render_rna
import svg
import sys

NODE_R = 6
PRIMARY_SPACE = 20
PAIR_SPACE = 23

CELL_PADDING = 40
TEXT_SIZE = 40

RENDER_IN_LETTERS = False

COLORS = {"r": [255, 0, 0],
          "g": [0, 155, 0],
          "b": [0, 0, 255],
          "k": [0, 0, 0],
          "y": [255, 255, 0],
          "c": [0, 255, 255],
          "m": [255, 0, 255],
          "w": [255, 255, 255],
          "grey": [100, 100, 100]}

def draw_rna(sequence, secstruct, colors, filename="secstruct"):
    r = render_rna.RNARenderer()

    pairmap = render_rna.get_pairmap_from_secstruct(secstruct)
    pairs = []
    for i in range(len(pairmap)):
        if pairmap[i] > i:
            pairs.append({"from":i, "to":pairmap[i], "p":1.0, "color":COLORS["grey"]})
    r.setup_tree(secstruct, NODE_R, PRIMARY_SPACE, PAIR_SPACE)
    size = r.get_size()

    cell_size = max(size) + CELL_PADDING * 2
    colors = [COLORS[x] for x in colors]

    svgobj = svg.svg("%s.svg" % filename, cell_size, cell_size)
    r.draw(svgobj, CELL_PADDING, CELL_PADDING, colors, pairs, sequence, RENDER_IN_LETTERS)

def main():
    with open(sys.argv[1]) as f:
        n = int(f.readline())
        for i in range(n):
            seq = f.readline().strip()
            secstruct = f.readline().strip()
            colorings = f.readline().strip().split(",")
            colors = []
            for coloring in colorings:
                [n, color] = coloring.split("x")
                colors += int(n) * [color]
            draw_rna(seq, secstruct, colors, "%s_%d" % (sys.argv[1], i))

if __name__ == "__main__":
    main()
