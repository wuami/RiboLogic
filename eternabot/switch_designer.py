import eterna_utils
import varna, draw_utils
import settings
import inv_utils
import random
import math
import ensemble_design
import unittest
import sys

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
    # ensure that secondary structures are the same length
    if(len(secstruct1) != len(secstruct2)):
        print "SS1 (%s) and SS2 (%s) lengths don't match" % (len(secstruct1), len(secstruct2))
        sys.exit(0)
    
    # generate pair mappings
    pairmap1 = eterna_utils.get_pairmap_from_secstruct(secstruct1)
    pairmap2 = eterna_utils.get_pairmap_from_secstruct(secstruct2)
    
    # +1 for each pair or single that doesn't match
    dist = 0
    udist = 0
    for ii in range(0,len(locks)):
        if(locks[ii] == "o"):
            continue
        elif(locks[ii] == "u"):
            if(secstruct1[ii] != secstruct2[ii]):
                udist += 1
        else:
            if(pairmap1[ii] != pairmap2[ii]):
                if(pairmap1[ii] > ii):
                    dist += 1
                if(pairmap2[ii] > ii):
                    dist += 1
    return dist + max(0, udist-threshold)

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

def rc(bases):
    rc = ""
    for base in reversed(bases):
        if base == "G":
            rc += "C"
        elif base == "C":
            rc += "G"
        elif base == "U":
            rc += "A"
        elif base == "A":
            rc += "U"
    return rc


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

class SwitchDesigner(object):

    def __init__(self, id, type, beginseq, constraints, targets, scoring_func, inputs = None, mode = "ghost"):
        # sequence information
        self.id = id
        self.type = type
        self.beginseq = beginseq
        self.sequence = beginseq
        self.n = len(self.sequence)
        self.constraints = constraints
        self.index_array = self.get_unconstrained_indices()

        # scoring
        self.scoring_func = scoring_func
        #if type == "multi_input_oligo":
        self.bp_distance_func = bp_distance_with_unpaired
        #else:
        #    self.bp_distance_func = bp_distance_with_constraint
        self.mode = mode

        # target information
        self.targets = targets
        self.n_targets = len(self.targets)
        self.inputs = inputs
        self.linker_length = 5
        self.linker = ensemble_design.get_sequence_string(["C" if i % 5 == 2 else "A" for i in range(self.linker_length)])
        #self.linker = "U"*self.linker_length
        self.design_linker = ""#GUUUCACCCCUAAACACCAC"
        self.oligotail = ""#AUUGUUAGUUAGGUAAAAAA"
        if type == "multi_input" or type == "multi_input_oligo":
            self.create_target_secstructs()

        self.update_sequence(*self.get_sequence_info(self.sequence))
        
        # maintain sequences
        self.update_best()
        self.all_solutions = []

    def create_target_secstructs(self):
        """ add oligo secstructs to target secstruct """
        for target in self.targets:
            if target['type'] == 'oligos':
                # create beginning linker
                target['secstruct'] += '.'*self.linker_length
                target['constrained'] += 'u'*self.linker_length
                if self.mode == "ghost":
                    target['fold_constraint'] = '.'*(self.n+self.linker_length)
                # add each input
                for i,input in enumerate(sorted(self.inputs)):
                    seq = self.inputs[input]
                    if input in target['inputs']:
                        if self.mode == "hairpin":
                            hairpin_len = (len(seq)-4)/2
                            secstruct = '('*hairpin_len + '.'*(len(seq)-2*hairpin_len) + ')'*hairpin_len
                        elif self.mode == "ghost":
                            secstruct = '.'*len(seq)
                            target['fold_constraint'] += 'x'*len(seq)
                        target['secstruct'] += secstruct
                        target['constrained'] += 'u'*len(seq)
                    else:
                        target['secstruct'] += '.'*len(seq)
                        target['constrained'] += 'o'*len(seq)
                        if self.mode == "ghost":
                            target['fold_constraint'] += '.'*len(seq)
                    # add linker sequence
                    #if i != len(self.inputs)-1:
                    target['secstruct'] += '.'*self.linker_length
                    target['constrained'] += 'u'*self.linker_length
                    if self.mode == "ghost":
                        target['fold_constraint'] += '.'*self.linker_length

    def get_unconstrained_indices(self):
        """
        get indices for positions that can change)
        """
        index_array = []
        for ii in range(0,self.n):
            if(self.constraints[ii] == "o"):
                index_array.append(ii)
        return index_array

    def get_solution(self):
        """
        return current best as solution
        """
        for i in range(self.n_targets):
            fold_sequence = self.get_fold_sequence(self.best_sequence, self.targets[i])
            print fold_sequence
            if self.mode == "ghost":
                print inv_utils.fold(fold_sequence, self.targets[i]['fold_constraint'])[0]
            else:
                print inv_utils.fold(fold_sequence)[0]
        return [self.best_sequence, self.best_bp_distance, self.best_design_score]

    def draw_solution(self):
        # initiate varna RNA visualizer
        v = varna.Varna()
        
        # get puzzle object and generate colormaps for each objective
        n = len(self.inputs)
        colormaps = draw_utils.get_colormaps(self.targets, self.inputs, self.n, self.linker_length, self.design_linker, n)
        
        # draw image for each sequence
        for i, target in enumerate(self.targets):
            filename = "%s/images/%s_%s-%s.png" % (settings.PUZZLE_DIR, self.id, 0, i)
            draw_utils.draw_secstruct_state(v, target, self.get_fold_sequence(self.sequence, target), colormaps[i], filename)

    def get_solutions(self):
        """
        return all possible solutions found
        """
        return self.all_solutions

    def get_random_solution(self):
        """
        return a random solution
        """
        r = random.randint(0, len(self.all_solutions)-1)
        solution = self.all_solutions[r]
        return [solution[0], 0, solution[1]]

    def get_design_score(self, sequence, secstruct):
        """
        calculates overall design score, which is sum of bp distance component and scoring function component
        """
        # in a small number of cases, get_design function causes error
        # if this happens, assume 0 score
        try:
            design = eterna_utils.get_design_from_sequence(sequence, secstruct[self.single_index])
            score = self.scoring_func(design)
            return score
        except:
            return 0

    def get_hairpin(self, seq, begin):
        """ turn sequence into hairpin """
        n = (len(seq)-4)/2
        if begin:
            return seq[0:n] + seq[n:len(seq)-n] + rc(seq[0:n])
        else:
            return rc(seq[-n:]) + seq[n:len(seq)-n] + seq[-n:]

    def get_fold_sequence(self, sequence, objective):
        """ append oligo sequences separated by & for type oligo """
        if objective['type'] == 'oligo':
            return '&'.join([sequence, objective['oligo_sequence']])
        elif objective['type'] == 'oligos':
            inputs = []
            for i,input in enumerate(sorted(self.inputs)):
                if input in objective['inputs']:
                    #inputs.append(self.inputs[input])
                    if self.mode == "ghost":
                        inputs.append(rc(self.inputs[input]))
                    elif i == len(self.inputs)-1:
                        inputs.append(self.get_hairpin(rc(self.inputs[input]), False))
                    else:
                        inputs.append(self.get_hairpin(rc(self.inputs[input]), True))
                else:
                    #inputs.append("U"*len(self.inputs[input]))
                    inputs.append(rc(self.inputs[input]))
            #input_sequence = self.linker.join(inputs)
            input_sequence = self.linker + self.linker.join(inputs) + self.linker
            return self.design_linker.join([sequence, input_sequence]) + self.oligotail
        else:
            return sequence 

    def get_sequence_info(self, sequence):
        """
        get sequence information - native fold, bp distance, pairmap, design
        for a particular sequence
        """
        native = []
        for i in range(self.n_targets):
            fold_sequence = self.get_fold_sequence(sequence, self.targets[i])
            if self.mode == "ghost":
                fold = inv_utils.fold(fold_sequence, self.targets[i]['fold_constraint'])[0]
            else:
                fold = inv_utils.fold(fold_sequence)[0]
            #if self.targets[i]['type'] == "oligos" and self.type != "multi_input_oligo":
            #    fold = fold[:len(sequence)]
            native.append(fold)
        bp_distance = self.score_secstructs(native)
        design_score = self.get_design_score(sequence, native)
        return [sequence, native, bp_distance, design_score]

    def reset_sequence(self):
        """
        reset sequence to the start sequence (for rerunning optimization)
        """
        self.sequence = self.beginseq
        self.update_sequence(*self.get_sequence_info(self.sequence))
        self.update_best()
        self.all_solutions = []

    def optimize_start_sequence(self):
        """
        optimize start sequence to include pairs, etc
        """
        secstruct = self.targets[self.single_index]['secstruct']
        constraints = convert_sequence_constraints(self.beginseq, self.constraints)
        sequence = ensemble_design.initial_sequence_with_gc_caps(secstruct, constraints, False)
        self.update_sequence(*self.get_sequence_info(sequence))
        return

    def update_sequence(self, sequence, native, bp_distance, score):
        """
        updates current sequence and related information
        """
        self.sequence = sequence
        [sequence, native, bp_distance, score] = self.get_sequence_info(sequence)
        self.native = native
        self.bp_distance = bp_distance
        self.design_score = score

    def update_best(self):
        """
        updates best to current sequence
        """
        self.best_sequence = self.sequence
        self.best_native = self.native
        self.best_bp_distance = self.bp_distance
        self.best_design_score = self.design_score

    def score_secstructs(self, secstruct):
        """
        calculates sum of bp distances for with and without oligo

        returns:
        sum of bp distances with and without the oligo 
        """
        distance = 0
        for i in range(self.n_targets):
            distance += self.bp_distance_func(secstruct[i], self.targets[i]['secstruct'], self.targets[i]['constrained'])
        return distance

    def check_secstructs(self, secstruct):
        """
        checks if current sequence matches target secondary structures
    
        return:
        boolean indicating whether the RNA folds to the targeted structure
        with and without the oligo
        """
        return self.score_secstructs(secstruct) == 0

    def check_current_secstructs(self):
        return self.score_secstructs(self.native) == 0

    def optimize_sequence(self, n_iterations, n_cool):
        """
        monte-carlo optimization of the sequence

        args:
        n_interations is the total number of iterations
        n_cool is the number of times to cool the system
        """
        bases = "GAUC"
        pairs = ["GC", "CG", "AU", "UA"]
    
        if len(self.index_array) == 0:
            return

        #print self.targets
        #self.optimize_start_sequence()

        T = 5

        def p_dist(dist, new_dist):
            """probability function"""
            if dist == 0:
               return 0
            return math.exp(-(new_dist-dist)/T)

        def p_score(score, new_score):
            """probability function for design scores"""
            return math.exp((new_score-score)/T)
    
        # loop as long as bp distance too large or design score too small
        for i in range(n_iterations):
            #random.shuffle(index_array)
            
            # pick random nucleotide in sequence
            rindex = self.index_array[int(random.random() * len(self.index_array))]
                
            mut_array = ensemble_design.get_sequence_array(self.sequence)
            mut_array[rindex] = ensemble_design.get_random_base()
            
            mut_sequence = ensemble_design.get_sequence_string(mut_array)
            [mut_sequence, native, bp_distance, score] = self.get_sequence_info(mut_sequence)

            # if current sequence is a solution, save to list
            if bp_distance == 0:
                self.all_solutions.append([mut_sequence, score])
            
            # if distance or score is better for mutant, update the current sequence
            if(random.random() < p_dist(self.bp_distance, bp_distance) or
               (bp_distance == self.bp_distance and random.random() < p_score(self.design_score, score))):
                self.update_sequence(mut_sequence, native, bp_distance, score)
                #print self.sequence, self.bp_distance
                #for i in range(self.n_targets):
                #    print self.native[i]
            
                # if distance or score is better for mutant than best, update the best sequence    
                if(bp_distance < self.best_bp_distance or
                   (bp_distance == self.bp_distance and score > self.best_design_score)):
                    self.update_best()

            # decrease temperature
            if i % n_iterations/n_cool == 0:
                T -= 0.1
                if T < 1:
                    T = 1
        
        return

class test_functions(unittest.TestCase):

    def setUp(self):
        self.puzzle = read_puzzle_json('switch_input.json')

    def test_check_secstructs(self):
        sequence = "ACAAGCUUUUUGCUCGUCUUAUACAUGGGUAAAAAAAAAACAUGAGGAUCACCCAUGUAAAAAAAAAAAAAAAAAAA"
        self.puzzle.update_sequence(*self.puzzle.get_sequence_info(sequence))
        self.assertTrue(self.puzzle.check_current_secstructs())

    def test_optimize_sequence(self):
        sequence = "AAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAACAUGAGGAUCACCCAUGUAAAAAAAAAAAAAAAAAAA"
        self.puzzle.update_sequence(*self.puzzle.get_sequence_info(sequence))
        self.puzzle.optimize_sequence(1000)
        print self.puzzle.get_solution()

