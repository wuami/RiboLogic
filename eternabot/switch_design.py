import inv_utils
import eterna_utils
import ensemble_design
import ensemble_utils
import unittest
import sys
import random
import json

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

class OligoPuzzle:

    def __init__(self, beginseq, constraints, target, scoring_func):
        # sequence information
        self.beginseq = beginseq
        self.sequence = beginseq
        self.n = len(self.sequence)
        self.constraints = constraints
        self.index_array = self.get_unconstrained_indices()

        self.scoring_func = scoring_func
        self.score_weight = 0.0001
        self.T = 100

        # target information
        self.target = target
        self.target_pairmap = []
        for x in target:
            self.target_pairmap.append(eterna_utils.get_pairmap_from_secstruct(x['secstruct']))
        self.n_targets = len(self.target)

        self.update_sequence(*self.get_sequence_info(self.sequence))
        
        # maintain best
        self.best_sequence = self.sequence
        self.best_native = self.native
        self.best_native_pairmap = self.native_pairmap
        self.best_design_score = self.design_score
        
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
        bp_distance = [] 
        for i in range(self.n_targets):
            bp_distance.append(bp_distance_with_constraint(self.target[i]['secstruct'],self.best_native[i],self.target[i]['constrained']))
        return [self.best_sequence, bp_distance, self.best_design_score]

    def get_design_score(self, sequence, secstruct):
        """
        calculates overall design score, which is sum of bp distance component and scoring function component
        """
        match_score = float(self.n - self.score_secstructs(secstruct))/self.n
        try:
            design = eterna_utils.get_design_from_sequence(sequence, secstruct[0])
            score = self.scoring_func(design)['finalscore']
            #score = 0
            #for i in range(self.n_targets):
            #    fold_sequence = self.get_fold_sequence(sequence, self.target[i])
            #    design = eterna_utils.get_design_from_sequence(fold_sequence, secstruct[i])
            #    score += self.scoring_func(design)['finalscore']
            return match_score + self.score_weight*score
        except IndexError:
            return match_score 

    def get_fold_sequence(self, sequence, objective):
        if objective['type'] == 'oligo':
            return '&'.join([sequence, objective['oligo_sequence']])
        else:
            return sequence 

    def get_sequence_info(self, sequence):
        """
        get sequence information - native fold, bp distance, pairmap, design
        for a particular sequence
        """
        native = []
        native_pairmap = []
        for x in self.target:
            fold_sequence = self.get_fold_sequence(sequence, x)
            native_current = inv_utils.fold(fold_sequence)[0]
            native.append(native_current)
            native_pairmap.append(eterna_utils.get_pairmap_from_secstruct(native_current))
        design_score = self.get_design_score(sequence, native)
        return [sequence, native, native_pairmap, design_score]

    def reset_sequence(self):
        self.sequence = self.beginseq
        self.update_sequence(*self.get_sequence_info(self.sequence))
        self.update_best()

    def update_sequence(self, sequence, native, native_pairmap, score):
        """
        updates current sequence and related information
        """
        self.sequence = sequence
        [sequence, native, native_pairmap, score] = self.get_sequence_info(sequence)
        self.native = native
        self.native_pairmap = native_pairmap
        self.design_score = score

    def update_best(self):
        """
        updates best to current sequence
        """
        self.best_sequence = self.sequence
        self.best_native = self.native
        self.best_native_pairmap = self.native_pairmap
        self.best_design_score = self.design_score

    def score_secstructs(self, secstruct):
        """
        calculates sum of bp distances for with and without oligo

        returns:
        sum of bp distances with and without the oligo 
        """
        distance = 0
        for i in range(self.n_targets):
            distance += bp_distance_with_constraint(secstruct[i], self.target[i]['secstruct'], self.target[i]['constrained'])
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

    def optimize_sequence(self, score_cutoff):
        """
        monte-carlo optimization of the sequence

        args:
        score_cutoff is the score at which to stop the loop
        """
        bases = "GAUC"
        pairs = ["GC", "CG", "AU", "UA"]
    
        if len(self.index_array) == 0:
            return
    
        # loop as long as bp distance too large or design score too small
        i = 0
        while(self.design_score < score_cutoff):
            i += 1
            #random.shuffle(index_array)
            
            # pick random nucleotide in sequence
            rindex = self.index_array[int(random.random() * len(self.index_array))]
                
            mut_array = ensemble_design.get_sequence_array(self.sequence)
            mut_array[rindex] = ensemble_design.get_random_base()
            
            mut_sequence = ensemble_design.get_sequence_string(mut_array)
            [mut_sequence, native, native_pairmap, score] = self.get_sequence_info(mut_sequence)
            
            # if distance or score is better for mutant, update the current sequence
            if(score > self.design_score or random.random() < score/self.design_score):
                self.update_sequence(mut_sequence, native, native_pairmap, score)
            
                # if distance or score is better for mutant than best, update the best sequence    
                if(score > self.best_design_score):
                    self.update_best()
            if i > 50000:
                print "optimization did not finish in 50000 iterations"
                break
        
        print "%s iterations" % i  
        return

def read_puzzle_json(filename):
    """
    read in puzzle as a json file
    """
    with open(filename, 'r') as f:
        p = json.loads(f.read())['data']['puzzle']
    beginseq = p['beginseq']
    constraints = p['locks']
    objective = json.loads(p['objective'])
    secstruct = [] 
    for o in objective:
        n = len(o['secstruct'])
        if 'structure_constrained_bases' not in o.keys() and 'anti_structure_constrained_bases' not in o.keys():
            constrained = ensemble_design.get_sequence_array('x'*n)
        else:
            constrained = ensemble_design.get_sequence_array('o'*n)
            struct = ensemble_design.get_sequence_array(o['secstruct'])
            if 'structure_constrained_bases' in o.keys() and len(o['structure_constrained_bases']) > 0:
                [lo, hi] = o['structure_constrained_bases']
                for i in range(lo, hi+1):
                    constrained[i] = 'x'
                del o['structure_constrained_bases']
            if 'anti_structure_constrained_bases' in o.keys() and len(o['anti_structure_constrained_bases']) > 0:
                [lo, hi] = o['anti_structure_constrained_bases']
                for i in range(lo, hi+1):
                    constrained[i] = 'x'
                    struct[i] = '.'
                del o['anti_secstruct'], o['anti_structure_constrained_bases']
            o['secstruct'] = ensemble_design.get_sequence_string(struct)
            o['constrained'] = ensemble_design.get_sequence_string(constrained)
        #secstruct[o['type']] = [struct, constrained]
        secstruct.append(o)
    strategy_names = ['example_gc60', 'penguian_clean_dotplot', 'berex_simplified_berex_test']
    ensemble = ensemble_utils.Ensemble("conventional", strategy_names, None)
    puzzle = OligoPuzzle(beginseq, constraints, secstruct, ensemble.score)
    return puzzle

def test_get_design(sequence):
    secstruct = inv_utils.fold(sequence)[0]
    try:
        eterna_utils.get_design_from_sequence(sequence, secstruct)
    except:
        print sequence, secstruct

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
        self.puzzle.optimize_sequence(1)
        print self.puzzle.get_solution()
        self.assertTrue(self.puzzle.check_current_secstructs())

def main():
    puzzle = read_puzzle_json("%s.json" % sys.argv[1])

    solutions = []
    scores = []
    n = 1
    i = 0
    while i < n:
        puzzle.reset_sequence()
        puzzle.optimize_sequence(1)
        assert puzzle.check_current_secstructs()
        sol = puzzle.get_solution()
        if sol[0] not in solutions:
            solutions.append(sol[0])
            scores.append(sol[2])
            print sol
            i += 1
        print "%s sequence calculated" % i

    with open(sys.argv[1] + ".out", 'w') as fout:
        for i in range(len(solutions)):
            fout.write("%s\t%1.6f\n" % (solutions[i], scores[i]))

if __name__ == "__main__":
    #unittest.main()
    main()


