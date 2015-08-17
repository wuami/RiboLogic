import sys

with open(sys.argv[1]) as f:
    lines = f.readlines()
    target = lines[5].split()[2]
    fold = lines[-2].split()[2]
    if target == fold:
        print lines[-1].split()[2]
