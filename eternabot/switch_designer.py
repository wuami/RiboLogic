import eterna_utils
import inv_utils
import random
import math
import ensemble_design
import unittest

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
        print "SS1 and SS2 lengths don't match"
        sys.exit(0)
    
    # generate pair mappings
    pairmap1 = eterna_utils.get_pairmap_from_secstruct(secstruct1)
    pairmap2 = eterna_utils.get_pairmap_from_secstruct(secstruct2)
    
    # +1 for each pair or single that doesn't match
    dist = 0
    
    for ii in range(0,len(pairmap1)):
        if(locks[ii] == "o"):
            continue
        if(pairmap1[ii] != pairmap2[ii]):
            if(pairmap1[ii] > ii):
                dist += 1
            if(pairmap2[ii] > ii):
                dist += 1
    
    return dist

def complement(base):
    if base == "G":
        return "C"
    elif base == "C":
        return "G"
    elif base == "U":
        return "A"
    elif base == "A":
        return "U"
    else:
        return "A"

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

class SwitchDesigner:

    def __init__(self, id, beginseq, constraints, targets, scoring_func, inputs = None):
        # sequence information
        self.id = id
        self.beginseq = beginseq
        self.sequence = beginseq
        self.n = len(self.sequence)
        self.constraints = constraints
        self.index_array = self.get_unconstrained_indices()

        self.scoring_func = scoring_func

        # target information
        self.targets = targets
        self.n_targets = len(self.targets)
        self.single_index = 0
        self.target_pairmap = []
        for i, target in enumerate(self.targets):
            if target['type'] == 'single':
                self.single_index = i
            self.target_pairmap.append(eterna_utils.get_pairmap_from_secstruct(target['secstruct']))
        self.inputs = inputs
        self.linker_length = 5
        self.linker = "U"*self.linker_length

        self.update_sequence(*self.get_sequence_info(self.sequence))
        
        # maintain sequences
        self.update_best()
        self.all_solutions = []
        
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
            print inv_utils.fold(fold_sequence)[0]
        return [self.best_sequence, self.best_bp_distance, self.best_design_score]

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
            score = self.scoring_func(design)['finalscore']
            return score
        except:
            return 0

    def get_fold_sequence(self, sequence, objective):
        """ append oligo sequences separated by & for type oligo """
        if objective['type'] == 'oligo':
            return '&'.join([sequence, objective['oligo_sequence']])
        elif objective['type'] == 'oligos':
            inputs = []
            for input in sorted(self.inputs):
                if input in objective['inputs']:
                    inputs.append(self.inputs[input])
                else:
                    inputs.append("U"*len(self.inputs[input]))
            input_sequence = self.linker.join(inputs)
            return '&'.join([sequence, input_sequence])
        else:
            return sequence 

    def get_sequence_info(self, sequence):
        """
        get sequence information - native fold, bp distance, pairmap, design
        for a particular sequence
        """
        native = []
        native_pairmap = []
        for i in range(self.n_targets):
            fold_sequence = self.get_fold_sequence(sequence, self.targets[i])
            fold = inv_utils.fold(fold_sequence)[0]
            if self.targets[i]['type'] == "oligos":
                fold = fold.split('&')[0]
            native.append(fold)
            native_pairmap.append(eterna_utils.get_pairmap_from_secstruct(native[i]))
        bp_distance = self.score_secstructs(native)
        design_score = self.get_design_score(sequence, native)
        return [sequence, native, native_pairmap, bp_distance, design_score]

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

    def update_sequence(self, sequence, native, native_pairmap, bp_distance, score):
        """
        updates current sequence and related information
        """
        self.sequence = sequence
        [sequence, native, native_pairmap, bp_distance, score] = self.get_sequence_info(sequence)
        self.native = native
        self.native_pairmap = native_pairmap
        self.bp_distance = bp_distance
        self.design_score = score

    def update_best(self):
        """
        updates best to current sequence
        """
        self.best_sequence = self.sequence
        self.best_native = self.native
        self.best_native_pairmap = self.native_pairmap
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
            distance += bp_distance_with_constraint(secstruct[i], self.targets[i]['secstruct'], self.targets[i]['constrained'])
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
            [mut_sequence, native, native_pairmap, bp_distance, score] = self.get_sequence_info(mut_sequence)

            # if current sequence is a solution, save to list
            if bp_distance == 0:
                self.all_solutions.append([mut_sequence, score])
            
            # if distance or score is better for mutant, update the current sequence
            if(random.random() < p_dist(self.bp_distance, bp_distance) or
               (bp_distance == self.bp_distance and random.random() < p_score(self.design_score, score))):
                self.update_sequence(mut_sequence, native, native_pairmap, bp_distance, score)
            
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

