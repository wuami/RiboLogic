import settings, fold_utils
import random, math
import re, string
import os, subprocess

def bp_distance_with_unpaired(secstruct1, secstruct2, locks, threshold=0):
    """
    calculates distance between two secondary structures
    
    args:
    secstruct1 is the first secondary structure
    secstruct2 is the second secondary structure
    locks specifies the positions that are constrained
    
    returns:
    bp distance between structures
    """
    order1 = ""
    if isinstance(secstruct1, list):
        order1 = secstruct1[1]
        secstruct1 = secstruct1[0]

    # ensure that secondary structures are the same length
    if(len(secstruct1) != len(secstruct2)):
        print "SS1 (%s) and SS2 (%s) lengths don't match" % (len(secstruct1), len(secstruct2))
        sys.exit(0)

    if not threshold:
        threshold = [[0,len(locks)-1,locks.count("u")]]
    
    # generate pair mappings
    if order1:
        pairmap1 = get_pairmap_from_secstruct([secstruct1, order1])
        ss_list = secstruct1.split('&')
        secstruct1 = "&".join([ss_list[order1.index(x)] for x in range(1,len(order1)+1)])
    else:
        pairmap1 = get_pairmap_from_secstruct(secstruct1)
    pairmap2 = get_pairmap_from_secstruct(secstruct2)
    
    # +1 for each pair or single that doesn't match
    dist = 0
    umatch = 0
    j = 0
    for i in range(0,len(locks)):
        if(locks[i] == "u"):
            if(secstruct1[i] == secstruct2[i]):
                umatch += 1
        elif locks[i] == "p":
            if secstruct1[i] in ['(', ')']:
                umatch += 1
        elif locks[i] == "x":
            if(pairmap1[i] != pairmap2[i]):
                if(pairmap1[i] > i):
                    dist += 1
                if(pairmap2[i] > i):
                    dist += 1
        else:
            continue
        if i == threshold[j][1]:
            dist += max(threshold[j][2]-umatch, 0)
            udist = 0
            if j != len(threshold)-1:
                j += 1
    return dist

class Scorer():
    def __init__(self, targets, nupack):
        self.MS2 = []
        self.indices = []
        self.nupack = nupack
        for target in targets:
            i = target['secstruct'].find('(((((.((....)))))))')
            if i == -1:
                self.MS2.append(False)
            else:
                self.MS2.append(True)
                self.indices = [i, i+18]
        if not any(self.MS2):
            self.scoring_func = self.pair_score
            for target in targets:
                pair_list = []
                stack = []
                for i in range(len(target['secstruct'])):
                    if target['constrained'][i] == "x":
                        if target['secstruct'][i] == "(":
                            stack.append(i)
                        elif target['secstruct'][i] == ")":
                            j = stack.pop()
                            pair_list.append([i,j])
                        elif target['secstruct'][i] == ".":
                            pair_list.append(-i)
                    elif target['constrained'][i] == "p":
                        pair_list.append(i)
                self.indices.append(pair_list)
        else:
            self.scoring_func = self.MS2_score

    def pair_score(self, sequences):
        score = 0.0
        for i, sequence in enumerate(sequences):
            n = len(sequence)
            dotplot = get_dotplot(sequence, self.nupack)
            pair_list = self.indices[i]
            for item in pair_list:
                if isinstance(item, list):
                    values = [pair[2] for pair in dotplot if (pair[0] == item[0] and pair[1] == item[1])]
                    if len(values) != 0:
                        score += sum(values)
                elif item < 0:
                    values = [pair[2] for pair in dotplot if (pair[0] == -item and pair[1] == n)]
                    if len(values) != 0:
                        score += sum(values)
                else:
                    values = [pair[2] for pair in dotplot if (item in pair and n not in pair)]
                    if len(values) != 0:
                        score += sum(values)
        return score

    def MS2_score(self, sequences):
        score = 0.0
        for i, sequence in enumerate(sequences):
            dotplot = get_dotplot(sequence, self.nupack)
            p = [pair for pair in dotplot if (pair[0] == self.indices[0] and pair[1] == self.indices[1])]
            if len(p) == 0:
                p = 0
            else:
                p = p[0][2]
            if self.MS2[i]:
                score += p
            else:
                score -= p
        return score*len(sequences[0])/len(sequences)

    def score(self, designs):
        return self.scoring_func(designs)
            
def get_bpp_scoring_func(targets, nupack): 
    s = Scorer(targets, nupack)
    return s.score


def get_dotplot(sequence, nupack=False, constraint=False):
    """ run ViennaRNA to get bp probability matrix """
    if nupack:
        return fold_utils.nupack_fold(sequence, nupack, bpp=True)
    filename = "".join(random.sample(string.lowercase,5))
    with open(filename+".fa",'w') as f:
        f.write(">%s\n" % filename)
        f.write("%s\n" % sequence)
        if constraint:
            options = " -C"
            f.write("%s\n" % constraint)
        else:
            options = ""
    if '&' in sequence:
        subprocess.call(os.path.join(settings.VIENNA_DIR,'RNAcofold') + options + ' -T 37.0 -p < %s' % filename + ".fa", shell=True, stdout=subprocess.PIPE)
    else:
        subprocess.call(os.path.join(settings.VIENNA_DIR,'RNAfold') + options + ' -T 37.0 -p < %s' % filename + ".fa", shell=True, stdout=subprocess.PIPE)

    # get info from output file
    try:
        file = open("%s_dp.ps" % filename, "r")
    except IOError:
        print "Can't find %s_dp.ps!" % filename
        sys.exit()

    dotps = file.read()
    file.close()
    os.remove(filename + ".fa")
    os.remove(filename + "_dp.ps")
    os.remove(filename + "_ss.ps")

    lines = re.findall('(\d+)\s+(\d+)\s+(\d*\.*\d*)\s+ubox',dotps)

    # create matrix containing index i, index j and pair probability
    dots = []

    for ii in range(0,len(lines)):
        dots.append([int(lines[ii][0]) - 1, int(lines[ii][1]) - 1, float(lines[ii][2])])

    return dots

def get_pairmap_from_secstruct(secstruct):
    """
    generates dictionary containing pair mappings

    args:
    secstruct contains secondary structure string

    returns:
    dictionary with pair mappings
    """
    order = ""
    if isinstance(secstruct, list):
        order = secstruct[1]
        secstruct = secstruct[0]
    pair_stack = []
    end_stack = []
    pk_pair_stack = []
    pk_end_stack = []
    pairs_array = []
    i_range = range(0,len(secstruct))

    # initialize all values to -1, meaning no pair
    for ii in i_range:
        pairs_array.append(-1)

    # assign pairs based on secstruct
    for ii in i_range:
        if(secstruct[ii] == "("):
            pair_stack.append(ii)
        elif(secstruct[ii] == ")"):
            if not pair_stack:
                end_stack.append(ii)
            else:
                index = pair_stack.pop()
                pairs_array[index] = ii
                pairs_array[ii] = index
        elif secstruct[ii] == "[":
            pk_pair_stack.append(ii)
        elif secstruct[ii] == "]":
            if not pk_pair_stack:
                pk_end_stack.append(ii)
            else:
                index = pk_pair_stack.pop()
                pairs_array[index] = ii
                pairs_array[ii] = index
    if len(pair_stack) == len(end_stack):
        n = len(pair_stack)
        for ii in range(n):
            pairs_array[pair_stack[ii]] = end_stack[-ii]
            pairs_array[end_stack[-ii]] = pair_stack[ii]
    else:
         print "ERROR: pairing incorrect %s" % secstruct

    # adjust pairs array for strand ordering
    if order:
        lengths = [len(x) for x in secstruct.split('&')]
        n = len(lengths)
        N = len(pairs_array)
        ordered_lengths = [lengths[order.index(i+1)] for i in range(n)]
        # shift indices
        new_pairs_array = pairs_array[:]
        for i in range(n):
            if order[i] != i+1:
                index_range = [sum(lengths[0:i]) + i, sum(lengths[0:i+1]) + i]
                actual_position = order[i]-1
                offset = sum(ordered_lengths[0:actual_position]) + actual_position
                new_pairs_array = [pairs_array[i]-index_range[0]+offset if pairs_array[i] >= index_range[0] and pairs_array[i] < index_range[1] else new_pairs_array[i] for i in range(N)]
        # reorder array
        pairs_array = []
        for i in range(n):
            pos = order.index(i+1)
            index_range = [sum(lengths[0:pos]) + pos, sum(lengths[0:pos+1]) + pos]
            pairs_array += new_pairs_array[index_range[0]:index_range[1]] + [-1]
        pairs_array = pairs_array[:-1]
            
    return pairs_array

def get_random_base(bases = "AUGC"): 
    """
    generates random base
    
    args:
    none
    
    returns:
    random base as char
    """
    nbases = len(bases)
    randn = int(math.floor(random.random() * nbases) % nbases)
    return bases[randn]

def get_rcs(base):
    if base == "G":
        return ["U", "C"]
    elif base == "U":
        return ["G", "A"]
    elif base == "A":
        return ["U"]
    elif base == "C":
        return ["G"]
    else:
        raise ValueError("invalid base: %s" % base)

def rc_single(base, pGU=0, possible_bases="AUGC"):
    complements = get_rcs(base) 
    if not any([b in possible_bases for b in complements]):
        raise ValueError("no complements to %s in %s" % (base, str(possible_bases)))
    if len(complements) > 1:
        if random.random() < pGU and complements[0] in possible_bases or complements[1] not in possible_bases:
            return complements[0]
        else:
            return complements[1]
    else:
        return complements[0] 

def rc(bases, pGU=0, possible_bases="AUGC"):
    rc = ""
    for base in bases:
        rc += rc_single(base, pGU, possible_bases)
    return rc[::-1]

def satisfies_constraints(sequence, beginseq, constraints):
    for i,letter in enumerate(constraints):
        if letter == "o":
            continue
        if beginseq[i] != sequence[i]:
            return False
    return True
