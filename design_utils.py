import settings, fold_utils
import random, math
import re, string
import os, sys, subprocess

def bp_distance(secstruct1, secstruct2, locks, threshold=0):
    """
    calculates distance between two secondary structures
    
    args:
    secstruct1 is the first secondary structure
    secstruct2 is the second secondary structure
    locks specifies the positions that are constrained
    
    returns:
    bp distance between structures
    """
    order = ''
    if isinstance(secstruct1, list):
        order = secstruct1[1]
        secstruct1 = secstruct1[0]

    # ensure that secondary structures are the same length
    if(len(secstruct1) != len(secstruct2)):
        print 'SS1 (%s) and SS2 (%s) lengths don\'t match' % (len(secstruct1), len(secstruct2))
        sys.exit(0)

    if not threshold:
        threshold = [[0,len(locks)-1,locks.count('u')]]
    
    # generate pair mappings
    if order:
        pairmap1 = get_pairmap_from_secstruct([secstruct1, order])
        ss_list = secstruct1.split('&')
        secstruct1 = '&'.join([ss_list[order.index(x)] for x in range(1,len(order)+1)])
    else:
        pairmap1 = get_pairmap_from_secstruct(secstruct1)
    pairmap2 = get_pairmap_from_secstruct(secstruct2)
    
    # +1 for each pair or single that doesn't match
    dist = 0
    match = 0
    j = 0
    nonmatching = []
    for i in range(0,len(locks)):
        if(locks[i] == 'u'):
            if(secstruct1[i] == secstruct2[i]):
                match += 1
        elif locks[i] == 'p':
            if secstruct1[i] in ['(', ')']:
                match += 1
        elif locks[i] == 'x':
            if(pairmap1[i] != pairmap2[i]):
                if(pairmap1[i] > i):
                    dist += 1
                    nonmatching.append(i)
                if(pairmap2[i] > i):
                    dist += 1
                    nonmatching.append(i)
        elif locks[i] == 'n':
            if(pairmap1[i] == pairmap2[i]):
                dist += 1
                nonmatching.append(i)
        else:
            continue
        if i == threshold[j][1]:
            dist += max(threshold[j][2]-match, 0)
            udist = 0
            if j != len(threshold)-1:
                j += 1
    return dist, nonmatching

class bpScorer():
    """
    class for base pair probability based scoring functions
    """
    def __init__(self, targets, nupack):
        self.nupack = nupack
        self.indices = []
        self.n = []

        # get list of base pairs to score
        for target in targets:
            self.n.append(len(target['secstruct']))
            pair_list = []
            unpair_list = []
            stack = []
            unstack = []
            for i in range(len(target['secstruct'])):
                if target['constrained'][i] == 'x':
                    if target['secstruct'][i] == '(':
                        stack.append(i)
                    elif target['secstruct'][i] == ')':
                        j = stack.pop()
                        pair_list.append([j,i])
                elif target['constrained'][i] == 'p':
                    pair_list.append(i)
                elif target['constrained'][i] == 'n':
                    if target['secstruct'][i] == '(':
                        unstack.append(i)
                    elif target['secstruct'][i] == ')':
                        j = unstack.pop()
                        unpair_list.append([j,i])
            self.indices.append([pair_list, unpair_list])

    def score(self, seq):
        dotplots = seq.bpps
        score = []
        for i, dotplot in enumerate(dotplots):
            n = self.n[i]
            pair_list, unpair_list = self.indices[i]
            for item in pair_list:
                if isinstance(item, list):
                    values = [pair[2] for pair in dotplot if set(pair[0:2]) == set(item[0:2])]
                    if len(values) != 0:
                        score.append(sum(values))
                else:
                    values = [pair[2] for pair in dotplot if (item in pair and n not in pair)]
                    if len(values) != 0:
                        score.append(sum(values))
            for item in unpair_list:
                values = [pair[2] for pair in dotplot if set(pair[0:2]) == set(item[0:2])]
                if len(values) != 0:
                    score.append(-sum(values))
        print score
        return sum(score)
            
def get_bpp_scoring_func(targets, nupack): 
    """
    get scorer for base pair probability based scoring function
    """
    s = bpScorer(targets, nupack)
    return s.score

def get_dotplot(sequence, nupack=False, constraint=False):
    """ run ViennaRNA to get bp probability matrix """
    if nupack:
        return fold_utils.nupack_fold(sequence, nupack, bpp=True)
    else:
        return fold_utils.vienna_fold(sequence, bpp=True)

def get_pairmap_from_secstruct(secstruct):
    """
    generates dictionary containing pair mappings

    args:
    secstruct contains secondary structure string

    returns:
    dictionary with pair mappings
    """
    order = ''
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
        if(secstruct[ii] == '('):
            pair_stack.append(ii)
        elif(secstruct[ii] == ')'):
            if not pair_stack:
                end_stack.append(ii)
            else:
                index = pair_stack.pop()
                pairs_array[index] = ii
                pairs_array[ii] = index
        elif secstruct[ii] == '[':
            pk_pair_stack.append(ii)
        elif secstruct[ii] == ']':
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
         print 'ERROR: pairing incorrect %s' % secstruct

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

def get_random_base(bases = 'AUGC'): 
    """
    generates random base
    """
    nbases = len(bases)
    randn = int(math.floor(random.random() * nbases) % nbases)
    return bases[randn]

def get_complements(base):
    """
    gets all possible complements of a base
    """
    if base == 'G':
        return ['U', 'C']
    elif base == 'U':
        return ['G', 'A']
    elif base == 'A':
        return ['U']
    elif base == 'C':
        return ['G']
    else:
        raise ValueError('invalid base: %s' % base)

def rc_single(base, pGU=0, possible_bases='AUGC'):
    """
    gets a complement of a single base, given possible choices and probability of GU pairs
    """
    complements = get_complements(base) 
    if not any([b in possible_bases for b in complements]):
        raise ValueError('no complements to %s in %s' % (base, str(possible_bases)))
    if len(complements) > 1:
        if random.random() < pGU and complements[0] in possible_bases or complements[1] not in possible_bases:
            return complements[0]
        else:
            return complements[1]
    else:
        return complements[0] 

def rc(bases, pGU=0, possible_bases='AUGC'):
    """
    gets reverse complement of a sequence
    """
    rc = ''
    for base in bases:
        rc += rc_single(base, pGU, possible_bases)
    return rc[::-1]

def satisfies_constraints(sequence, beginseq, constraints):
    """
    determines if given sequence satisfies sequence constraints
    """
    for i,letter in enumerate(constraints):
        if letter == 'o':
            continue
        if beginseq[i] != sequence[i]:
            return False
    return True

def weighted_choice(seq, w):
    r = random.uniform(0, sum(w))
    current = 0
    for i, choice in enumerate(seq):
        if current + w[i] >= r:
            return choice
        current += w[i]
