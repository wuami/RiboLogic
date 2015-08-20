import sys, os
import random, string
import vienna_parameters
import subprocess
import eterna_utils
import random
import re
from subprocess import Popen, PIPE, STDOUT, check_call
import thread, time, sys
from threading import Timer
import settings


DEFAULT_TEMPERATURE = 37.0
BASES = ['A','U','G','C']
#fold a sequence
#@param seq:sequence
#@return: [parenthesis notation, energy]

def fold(seq, cotransc=False, constraint=False):
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
        p = Popen([os.path.join(settings.VIENNA_DIR,'CoFold'), '--distAlpha', '0.5', '--distTau', '640', '--noPS', options], stdout=PIPE, stdin=PIPE, stderr=STDOUT)
    elif '&' in seq:
        p = Popen([os.path.join(settings.VIENNA_DIR,'RNAcofold'), '-T','37.0', '-noPS', options], stdout=PIPE, stdin=PIPE, stderr=STDOUT)
    else:
        p = Popen([os.path.join(settings.VIENNA_DIR,'RNAfold'), '-T','37.0', '-noPS', options], stdout=PIPE, stdin=PIPE, stderr=STDOUT)
    pair= p.communicate(input=input)[0]
    p.wait()

    # split result by whitespace
    toks = re.split('\s+| \(?\s?',pair)
    ret= []
    ret.append(toks[1])
    ret.append(toks[2][1:-1])
    return ret

def nupack_fold(seq, bpp = False):
    """
    folds sequence using nupack
    """
    rand_string = ''.join(random.choice(string.ascii_lowercase) for _ in range(5))
    split = seq.split("&")
    f = open("%s.in" % rand_string, "w")
    f.write("%s\n" % len(split))
    for seq in split:
        f.write("%s\n" % seq)
    f.write("1\n")
    os.system("cp %s/%s.list %s.list" % (settings.NUPACK_DIR, len(split), rand_string))
    f.close()
    options = ['-material', 'rna1999', '-ordered', '-mfe']#, '-quiet']
    if bpp:
        options.append('-pairs')
    p = Popen([os.path.join(settings.NUPACK_DIR,'complexes')] + options + [rand_string], stdout=PIPE, stdin=PIPE, stderr=STDOUT)
    p.wait()
    # get strand ordering
    strands_dict = {}
    with open("%s.ocx-key" % rand_string) as f_key:
        for line in f_key:
            if not line.startswith("%"):
                spl = line.strip().split("\t")
                strands_dict[(spl[0], spl[1])] = spl[2:]
    # get mfe structure
    with open("%s.ocx-mfe" % rand_string) as f_mfe:
        best_complex = -1
        best_strands = []
        best_mfe = float("inf")
        best_fold = ""
        line = f_mfe.readline()
        while line:
            if line.startswith("% complex"):
                m = re.search(r"complex([0-9]+)-order([0-9]+)", line)
                complex = (m.group(1), m.group(2))
                strands = strands_dict[complex]
                if len(strands) == len(set(strands)) and len(strands) == len(split):
                    f_mfe.readline()
                    energy = float(f_mfe.readline().strip())
                    if energy < best_mfe:
                        best_complex = complex
                        best_strands = [int(x) for x in strands]
                        best_mfe = energy
                        best_fold = f_mfe.readline().strip()
            line = f_mfe.readline()

    if bpp:
        bpp_matrix = []
        with open("%s.cx-epairs" % rand_string) as f_pairs:
            line = f_pairs.readline()
            while line:
                if line.startswith("% complex%s" % complex[0]):
                    f_pairs.readline()
                    line = f_pairs.readline()
                    while not line.startswith("%"):
                        bpp_matrix.append(line.strip().split())
        os.system("rm %s*" % rand_string)
        return bpp_matrix

    # get secondary structures in order
    os.system("rm %s*" % rand_string)
    return [best_fold.replace("+", "&"), best_mfe, best_strands]
    secstructs = best_fold.split("+")
    folds = []
    for i in range(1, len(split)+1):
        if i in best_strands:
            folds.append(secstructs[best_strands.index(i)])
        else:
            folds.append("."*len(split[i-1]))
    return ["&".join(folds), best_mfe]

def fill_gc(elem , pair_map , seq, rand ):
    if(elem.type_ != eterna_utils.RNAELEMENT_STACK):
        return
    indices = elem.indices_
    length = len(indices)
    for ii in range(0,length):
        idx = indices[ii]
        if(pair_map[idx]<idx):
            continue
        if(rand.randint(0,1)==0):
            seq[idx]="G"
            seq[pair_map[idx]]="C"
        else:
            seq[idx]="C"
            seq[pair_map[idx]]="G"

def tout():
    thread.interrupt_main()

def timeout(func, args=(), timeout_duration=10, default=[]):
    ret = []
    finished=False
    seq=default
    t0=timeout
    try:
        timer = Timer(timeout_duration, tout)
        timer.start()
        t0=time.clock()
        seq = func(*args)
        finished=True
        t0=time.clock()-t0
        timer.cancel()
    except:
        print "time out"
    if(finished):
        ret.append(seq)
        ret.append(t0)
    else:
        ret.append(default)
        ret.append("timeout")
    return ret

