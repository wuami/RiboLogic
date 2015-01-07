import inv_utils
import eterna_utils
import ensemble_design
import ensemble_utils
import unittest
import sys
import random
import math

def bp_distance_with_constraint(secstruct1, secstruct2, constraint):
    """
    calculates distance between two secondary structures
    
    args:
    secstruct1 is the first secondary structure
    secstruct2 is the second secondary structure
    constraint specifies the positions that are constrained
    
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
        if(constraint[ii] == "N"):
            continue
        if(pairmap1[ii] != pairmap2[ii]):
            if(pairmap1[ii] > ii):
                dist += 1
            if(pairmap2[ii] > ii):
                dist += 1
    
    return dist

class OligoPuzzle:

    def __init__(self, constraints, oligo, target, scoring_func):
        # sequence information
        self.sequence = constraints.replace("N", "A")
        self.n = len(self.sequence)
        self.constraints = constraints
        self.oligo = oligo

        self.scoring_func = scoring_func
        self.score_weight = 0
        self.T = 100

        # target information
        self.target = target
        self.target_pairmap = {}
        self.target_pairmap['oligo']  = eterna_utils.get_pairmap_from_secstruct(self.target['oligo'][0])
        self.target_pairmap['no_oligo']  = eterna_utils.get_pairmap_from_secstruct(self.target['no_oligo'][0])

        self.update_sequence(self.sequence)
        
        # maintain best
        self.best_sequence = self.sequence
        self.best_native = self.native
        self.best_native_pairmap = self.native_pairmap
        #self.best_design = self.design
        self.best_design_score = self.design_score
        
        self.index_array = self.get_unconstrained_indices()

    def get_unconstrained_indices(self):
        """
        get indices for positions that can change
        """
        index_array = []
        for ii in range(0,self.n):
            if(self.constraints[ii] == "N"):
                index_array.append(ii)
        return index_array

    def get_solution(self):
        """
        return current best as solution
        """
        bp_distance = {}
        bp_distance['no_oligo'] = bp_distance_with_constraint(self.target['no_oligo'][0],self.best_native['no_oligo'],self.target['no_oligo'][1])
        bp_distance['oligo'] = bp_distance_with_constraint(self.target['oligo'][0],self.best_native['oligo'],self.target['oligo'][1])
        return [self.best_sequence, bp_distance, self.best_design_score]

    def get_design_score(self, secstruct):#, design):
        """
        calculates overall design score, which is sum of bp distance component and scoring function component
        """
        return (self.n - self.score_secstructs(secstruct))# + self.score_weight*self.scoring_func(design)['finalscore']

    def get_sequence_info(self, sequence):
        """
        get sequence information - native fold, bp distance, pairmap, design
        for a particular sequence
        """
        native = {}
        native['no_oligo'] = inv_utils.fold(sequence)[0]
        native['oligo'] = inv_utils.cofold('&'.join([sequence, self.oligo]))[0][:self.n]
        native_pairmap = {}
        native_pairmap['oligo'] = eterna_utils.get_pairmap_from_secstruct(native['oligo'])
        native_pairmap['no_oligo'] = eterna_utils.get_pairmap_from_secstruct(native['no_oligo'])
        #design = eterna_utils.get_design_from_sequence(sequence,native['no_oligo'])
        design_score = self.get_design_score(native)#, design)
        return [native, native_pairmap, design_score]#, design]

    def reinitialize_sequence(self):
        sequence = self.constraints.replace("N", "A")
        self.update_sequence(sequence)

    def update_sequence(self, sequence):
        """
        updates current sequence and related information
        """
        self.sequence = sequence
        [native, native_pairmap, score] = self.get_sequence_info(sequence)
        self.native = native
        self.native_pairmap = native_pairmap
        #self.design = design
        self.design_score = score

    def update_best(self):
        """
        updates best to current sequence
        """
        self.best_sequence = self.sequence
        self.best_native = self.native
        self.best_native_pairmap = self.native_pairmap
        #self.best_design = self.design
        self.best_design_score = self.design_score

    def score_secstructs(self, secstruct):
        """
        calculates sum of bp distances for with and without oligo

        returns:
        sum of bp distances with and without the oligo 
        """
        return bp_distance_with_constraint(secstruct['oligo'], self.target['oligo'][0], self.target['oligo'][1]) + \
               bp_distance_with_constraint(secstruct['no_oligo'], self.target['no_oligo'][0], self.target['no_oligo'][1])

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
    
        # initialize variables for iteration
        walk_iter = 0
        stale_move = 0
    
        if len(self.index_array) == 0:
            return
    
        # loop as long as bp distance too large or design score too small
        while(self.design_score < score_cutoff):
            #random.shuffle(index_array)
            
            # pick random nucleotide in sequence
            rindex = self.index_array[int(random.random() * len(self.index_array))]
                
            mut_array = ensemble_design.get_sequence_array(self.sequence)
            mut_array[rindex] = ensemble_design.get_random_base()
            
            mut_sequence = ensemble_design.get_sequence_string(mut_array)
            [native, native_pairmap, score] = self.get_sequence_info(mut_sequence)
            
            # if distance or score is better for mutant, update the current sequence
            if(score > self.design_score or random.random() < score/self.design_score):
                self.update_sequence(mut_sequence)
            
                # if distance or score is better for mutant than best, update the best sequence    
                if(score > self.best_design_score):
                    self.update_best()
                    
        return

class test_functions(unittest.TestCase):

    def setUp(self):
        with open('switch_input.txt', 'r') as f:
            secstruct = {}
            secstruct['oligo'] = f.readline().split()
            secstruct['no_oligo'] = f.readline().split()
            constraints = f.readline().strip()
            oligo_sequence = f.readline().strip()
        strategy_names = ['example_gc60', 'penguian_clean_dotplot', 'berex_simplified_berex_test']
        ensemble = ensemble_utils.Ensemble("conventional", strategy_names, None)
        self.puzzle = OligoPuzzle(constraints, oligo_sequence, secstruct, ensemble.score)

    def test_check_secstructs(self):
        sequence = "ACAAGCUUUUUGCUCGUCUUAUACAUGGGUAAAAAAAAAACAUGAGGAUCACCCAUGUAAAAAAAAAAAAAAAAAAA"
        self.puzzle.update_sequence(sequence)
        self.assertTrue(self.puzzle.check_current_secstructs())

    def test_optimize_sequence(self):
        self.puzzle.update_sequence("AAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAACAUGAGGAUCACCCAUGUAAAAAAAAAAAAAAAAAAA")
        self.puzzle.optimize_sequence(77)
        print self.puzzle.get_solution()
        self.assertTrue(self.puzzle.check_current_secstructs())

def main():
    with open('switch_input.txt', 'r') as f:
        secstruct = {}
        secstruct['oligo'] = f.readline().split()
        secstruct['no_oligo'] = f.readline().split()
        constraints = f.readline().strip()
        oligo_sequence = f.readline().strip()
    strategy_names = ['example_gc60', 'penguian_clean_dotplot', 'berex_simplified_berex_test']
    ensemble = ensemble_utils.Ensemble("conventional", strategy_names, None)
    puzzle = OligoPuzzle(constraints, oligo_sequence, secstruct, ensemble.score)

    fout = open('switch_output.txt', 'w')
    for i in xrange(100):
        puzzle.reinitialize_sequence()
        puzzle.optimize_sequence(77)
        assert puzzle.check_current_secstructs()
        fout.write("%s\n" % puzzle.get_solution()[0])
        print "%s sequence calculated" % i
    fout.close()

if __name__ == "__main__":
    #unittest.main()
    main()


