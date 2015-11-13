import subprocess, re, string, random, os, settings

def vienna_fold(seq, cotransc=False, constraint=False):
    """
    folds sequence using Vienna

    args:
    seq is the sequence string

    returns:
    secondary structure
    """
    # run ViennaRNA
    if constraint:
        options = "-C"
        input = seq + "\n" + constraint
    else:
        options = ""
        input = ''.join(seq)
    if cotransc:
        p = subprocess.Popen([os.path.join(settings.VIENNA_DIR,'CoFold'), '--distAlpha', '0.5', '--distTau', '640', '--noPS', options], stdout=subprocess.PIPE, stdin=subprocess.PIPE, stderr=subprocess.STDOUT)
    elif '&' in seq:
        p = subprocess.Popen([os.path.join(settings.VIENNA_DIR,'RNAcofold'), '-T','37.0', '-noPS', options], stdout=subprocess.PIPE, stdin=subprocess.PIPE, stderr=subprocess.STDOUT)
    else:
        p = subprocess.Popen([os.path.join(settings.VIENNA_DIR,'RNAfold'), '-T','37.0', '-noPS', options], stdout=subprocess.PIPE, stdin=subprocess.PIPE, stderr=subprocess.STDOUT)
    pair= p.communicate(input=input)[0]
    p.wait()

    # parse the result
    toks = re.search('([AUGC]+)\s*([\.\)\(]+)\s+\(\s*([-0-9\.]+)\s*\)', pair)
    return [toks.group(2), float(toks.group(3))]

def nupack_fold(seq, oligo_conc, bpp = False):
    """
    finds most prevalent structure using nupack partition function
    """
    os.system("source ~/.bashrc")
    rand_string = ''.join(random.choice(string.ascii_lowercase) for _ in range(6))
    split = seq.split("&")
    with open("%s.in" % rand_string, "w") as f:
        f.write("%s\n" % len(split))
        for seq in split:
            f.write("%s\n" % seq)
        f.write("1\n")
    os.system("cp %s/%s.list %s.list" % (settings.NUPACK_DIR, len(split), rand_string))
    with open("%s.con" % rand_string, "w") as f:
        if isinstance(oligo_conc, list):
            f.write("%s\n" % "\n".join([str(x) for x in oligo_conc]))
        else:
            f.write("%s\n" % oligo_conc * (len(split)-1))
        f.write("1e-9\n")
    options = ['-material', 'rna', '-ordered', '-mfe']#, '-quiet']
    if bpp:
        options.append('-pairs')
    p = subprocess.Popen([os.path.join(settings.NUPACK_DIR,'complexes')] + options + [rand_string], stdout=subprocess.PIPE, stdin=subprocess.PIPE, stderr=subprocess.PIPE)
    p.wait()
    p = subprocess.Popen([os.path.join(settings.NUPACK_DIR,'concentrations'), '-ordered', rand_string], stdout=subprocess.PIPE, stdin=subprocess.PIPE, stderr=subprocess.PIPE)
    p.wait()
    # get mfe
    with open("%s.eq" % rand_string) as f_eq:
        for line in f_eq:
            if not line.startswith("%"):
                complex = line.strip().split()
                if int(complex[len(split)+1]) and not bpp:
                    break
                if all([int(x) for x in complex[2:(2+len(split))]]) and bpp:
                    break
    # get strand ordering
    with open("%s.ocx-key" % rand_string) as f_key:
        for line in f_key:
            if line.startswith("%s\t%s" % (complex[0], complex[1])):
                strands = [int(x) for x in line.strip().split()[2:]]
                break
    # get mfe structure
    with open("%s.ocx-mfe" % rand_string) as f_mfe:
        line = f_mfe.readline()
        while line:
            if line.startswith("%% complex%s-order%s" % (complex[0], complex[1])):
                f_mfe.readline()
                energy = f_mfe.readline().strip()
                secstruct = f_mfe.readline().strip()
                break
            line = f_mfe.readline()

    if bpp:
        bpp_matrix = []
        with open("%s.cx-epairs" % rand_string) as f_pairs:
            line = f_pairs.readline()
            while line:
                if line.startswith("%% complex%s" % complex[0]):
                    f_pairs.readline()
                    line = f_pairs.readline()
                    while not line.startswith("%"):
                        bp = line.strip().split()
                        bpp_matrix.append([int(bp[0]), int(bp[1]), float(bp[2])])
                        line = f_pairs.readline()
                    break
                line = f_pairs.readline()
        os.system("rm %s*" % rand_string)
        return bpp_matrix

    # get full secondary structure
    for i in range(len(split)):
        if i+1 not in strands:
            secstruct += "&" + "."*len(split[i])
            strands.append(i+1)
    os.system("rm %s*" % rand_string)
    return [secstruct.replace("+", "&"), energy, strands]
