import eterna_utils
import sys
import math, random
import settings, os

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
        pairmap1 = eterna_utils.get_pairmap_from_secstruct([secstruct1, order1])
        ss_list = secstruct1.split('&')
        secstruct1 = "&".join([ss_list[order1.index(x)] for x in range(1,len(order1)+1)])
    else:
        pairmap1 = eterna_utils.get_pairmap_from_secstruct(secstruct1)
    pairmap2 = eterna_utils.get_pairmap_from_secstruct(secstruct2)
    
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
                dist += 1
                #if(pairmap1[i] > i):
                #    dist += 1
                #if(pairmap2[i] > i):
                #    dist += 1
        else:
            continue
        if i == threshold[j][1]:
            dist += max(threshold[j][2]-umatch, 0)
            udist = 0
            if j != len(threshold)-1:
                j += 1
    return dist

def bp_distance_with_constraint(secstruct1, secstruct2, locks):
    """
    calculates distance between two secondary structures
    
    args:
    secstruct1 is the first secondary structure
    secstruct2 is the second secondary structure
    locks specifies the positions that are constrained
    
    returns:
    bp distance between structures
    """
    # ensure that secondary structures are the same length
    if(len(secstruct1) != len(secstruct2)):
        print "SS1 (%s) and SS2 (%s) lengths don't match" % (len(secstruct1), len(secstruct2))
        sys.exit(0)
    
    # generate pair mappings
    pairmap1 = eterna_utils.get_pairmap_from_secstruct(secstruct1)
    pairmap2 = eterna_utils.get_pairmap_from_secstruct(secstruct2)
    
    # +1 for each pair or single that doesn't match
    dist = 0
    for ii in range(0,len(locks)):
        if(locks[ii] == "o"):
            continue
        if(pairmap1[ii] != pairmap2[ii]):
            if(pairmap1[ii] > ii):
                dist += 1
            if(pairmap2[ii] > ii):
                dist += 1
    return dist

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

def rc_single(base, pGU=0.2, possible_bases="AUGC"):
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

def rc(bases, pGU=0.2, possible_bases="AUGC"):
    rc = ""
    for base in bases:
        rc += rc_single(base, pGU, possible_bases)
    return rc

def get_different_base(base):
    bases ="GCAU"
    randn = int(math.floor(random.random() * 4) % 4)
    while bases[randn] == base or (base == "C" and bases[randn] == "U"):
        randn = int(math.floor(random.random() * 4) % 4)
    return bases[randn]

def satisfies_constraints(sequence, beginseq, constraints):
    for i,letter in enumerate(constraints):
        if letter == "o":
            continue
        if beginseq[i] != sequence[i]:
            return False
    return True

def convert_sequence_constraints(sequence, constraints):
    """
    convert "xo" constraints to "N" constraint format
    """
    result = ""
    for i,letter in enumerate(constraints):
        if letter == 'o':
            result += 'N'
        else:
            result += sequence[i]
    return result

def get_ensemble_scoring_func():
    """ create eternabot ensemble scoring function """
    strategy_names = ['merryskies_only_as_in_the_loops', 'aldo_repetition', 'dejerpha_basic_test', 'eli_blue_line', 'clollin_gs_in_place', 'quasispecies_test_by_region_boundaries', 'eli_gc_pairs_in_junction', 'eli_no_blue_nucleotides_in_hook', 'mat747_31_loops', 'merryskies_1_1_loop', 'xmbrst_clear_plot_stack_caps_and_safe_gc', 'jerryp70_jp_stratmark', 'eli_energy_limit_in_tetraloops', 'eli_double_AUPair_strategy', 'eli_green_blue_strong_middle_half', 'eli_loop_pattern_for_small_multiloops', 'eli_tetraloop_similarity', 'example_gc60', 'penguian_clean_dotplot', 'eli_twisted_basepairs', 'aldo_loops_and_stacks', 'eli_direction_of_gc_pairs_in_multiloops_neckarea', 'eli_multiloop_similarity', 'eli_green_line', 'ding_quad_energy', 'quasispecies_test_by_region_loops', 'berex_berex_loop_basic', 'eli_legal_placement_of_GUpairs', 'merryskies_1_1_loop_energy', 'ding_tetraloop_pattern', 'aldo_mismatch', 'eli_tetraloop_blues', 'eli_red_line', 'eli_wrong_direction_of_gc_pairs_in_multiloops', 'deivad_deivad_strategy', 'eli_direction_of_gc_pairs_in_multiloops', 'eli_no_blue_nucleotides_strategy', 'berex_basic_test', 'eli_numbers_of_yellow_nucleotides_pr_length_of_string', 'kkohli_test_by_kkohli']
    weights_file_name = "no_validation_training/weights_sparse_5.overall.txt"
    scores_file_name = "no_validation_training/predicted_score_sparse_5.overall.unnormalized.txt"
    weights_f = open(os.path.join(settings.RESOURCE_DIR, weights_file_name),"r")
    weights = []
    for line in weights_f:
        weights.append(float(line))
    ensemble = ensemble_utils.Ensemble("sparse", strategy_names, weights)
    def scoring_func(designs):
        return ensemble.score(designs[0])['finalscore']
    return scoring_func

class Scorer():
    def __init__(self, targets):
        self.MS2 = []
        self.indices = []
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

    def pair_score(self, designs):
        score = 0.0
        for i, design in enumerate(designs):
            n = len(design['sequence'])
            pair_list = self.indices[i]
            for item in pair_list:
                if isinstance(item, list):
                    values = [pair[2] for pair in design['dotplot'] if (pair[0] == item[0] and pair[1] == item[1])]
                    if len(values) != 0:
                        score += sum(values)
                elif item < 0:
                    values = [pair[2] for pair in design['dotplot'] if (pair[0] == -item and pair[1] == n)]
                    if len(values) != 0:
                        score += sum(values)
                else:
                    values = [pair[2] for pair in design['dotplot'] if (item in pair and n not in pair)]
                    if len(values) != 0:
                        score += sum(values)
        return score

    def MS2_score(self, designs):
        score = 0.0
        for i,design in enumerate(designs):
            p = [pair for pair in design['dotplot'] if (pair[0] == self.indices[0] and pair[1] == self.indices[1])]
            if len(p) == 0:
                p = 0
            else:
                p = p[0][2]
            if self.MS2[i]:
                score += p
            else:
                score -= p
        return score*len(designs[0]['sequence'])/len(designs)

    def score(self, designs):
        return self.scoring_func(designs)
            
def get_bpp_scoring_func(targets): 
    s = Scorer(targets)
    return s.score

def get_strategy_scoring_func(name):
    strategy = eterna_utils.load_strategy_from_file(os.path.join(settings.STRATEGY_DIR, name + '.py'))
    return strategy.score
