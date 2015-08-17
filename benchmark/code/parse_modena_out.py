import sys

max = -1
bestseq = None
with open(sys.argv[1]) as f:
    line = f.readline()
    while line:
        line = line.strip()
        if line.startswith("Individual"):
            spl = line.split()
            if float(spl[-1]) == 1.0:
                eps = float(spl[-2])
                if eps > max:
                    max = eps
                    f.readline()
                    bestseq = f.readline().strip()
        line = f.readline()

if max > 0:
    print bestseq
