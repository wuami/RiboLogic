import design_utils, fold_utils
import sequence_graph
import varna, draw_utils
import settings
import math, random
import unittest
import sys
import multiprocessing

class SwitchDesigner(object):

    def __init__(self, id, type, beginseq, constraints, targets, **kwargs):
        # sequence information
        self.id = id
        self.type = type
        self.beginseq = beginseq
        self.sequence = beginseq
        self.n = len(self.sequence)
        self.constraints = constraints

        self.mode = kwargs.get("mode", "nupack")
        add_rcs = kwargs.get("add_rcs", False)
        self.strandbonus = kwargs.get("strandbonus", False)
        self.print_ = kwargs.get("print_", False)
        self.inputs = kwargs.get("inputs", {})
        self.aptamer = False

        self.targets = self.parse_targets(targets)
        self.sequence_graph = sequence_graph.SequenceGraph(self.inputs, targets, constraints, beginseq, add_rcs, False, autocomplement=True)
        self.target_oligo_conc = 1e-7

        # scoring
        self.scoring_func = design_utils.get_bpp_scoring_func(targets, self.mode == "nupack")

        #if type == "multi_input_oligo":
        self.bp_distance_func = design_utils.bp_distance
        self.greedy = False
        self.oligo_conc = 1e7

        # target information
        self.n_targets = len(self.targets)
        n_cores = min(16, self.n_targets)
        
        # print puzzle info
        print self.constraints
        for i, target in enumerate(targets):
            print "-> state %d" % i
            print self.get_fold_sequence(self.sequence, target)
            print target['secstruct']
            print target['constrained']
            if 'fold_constraint' in target:
                print target['fold_constraint']

        self.update_current(*self.get_sequence_info(self.sequence))
        
        # maintain sequences
        self.update_best()
        self.all_solutions = []

    def parse_targets(self, targets):
        """
        generate full secondary structure and constraint strings
        """
        for target in targets:
            secstruct = ""
            constrained = ""
            if target['type'] == 'oligos' and len(target['inputs']) > 0 and '&' not in target['secstruct']:
                for input in target['inputs']:
                    secstruct += '.'*len(self.inputs[input]) + '&'
                    constrained += 'o'*len(self.inputs[input]) + 'x'
                target['secstruct'] = secstruct + target['secstruct']
                target['constrained'] = constrained + target['constrained']
            elif target['type'] == 'aptamer':
                self.aptamer = float(target['concentration'])
                fold_constraint = list(target['secstruct'])
                for i, fold in enumerate(fold_constraint):
                    if i in target['site'] and fold == ".":
                        fold_constraint[i] = "x"
                    elif i not in target['site']:
                        fold_constraint[i] = "."
                target['fold_constraint'] = "".join(fold_constraint)
            if '&' in target['secstruct']:
                breaks = [i for i, char in enumerate(target['secstruct']) if char == '&']
                constrained = list(target['constrained'])
                for i in breaks:
                    constrained[i] = 'x'
                target['constrained'] = "".join(constrained)
        return targets

    def get_fold_sequence(self, sequence, target):
        """ append oligo sequences separated by & for type oligo """
        if target['type'] == 'oligos':
            return '&'.join([self.inputs[x] for x in sorted(target['inputs'])] + [sequence])
        else:
            return sequence 

    def get_solution(self):
        """
        return current best as solution
        """
        for i in range(self.n_targets):
            fold_sequence = self.get_fold_sequence(self.best_sequence, self.targets[i])
            print fold_sequence
            if self.mode == "vienna":
                if self.targets[i]['type'] == "aptamer":
                    print fold_utils.vienna_fold(fold_sequence, self.targets[i]['fold_constraint'])[0]
                else:
                    print fold_utils.vienna_fold(fold_sequence)[0]
            else:
                print fold_utils.nupack_fold(fold_sequence, self.target_oligo_conc*self.oligo_conc)[0]
        return [self.best_sequence, self.best_bp_distance, self.best_design_score]

    def draw_solution(self, name):
        """ draw each state for current solution """
        # initiate varna RNA visualizer
        v = varna.Varna()
        
        # get puzzle object and generate colormaps for each target
        n = len(self.inputs)
        colormap = draw_utils.get_colormaps(self.targets, self.inputs, self.n, self.linker_length, self.design_linker, n)
        
        # draw image for each condition
        for i, target in enumerate(self.targets):
            filename = "%s/images/%s_%s-%s.png" % (settings.PUZZLE_DIR, self.id, name, i)
            draw_utils.draw_secstruct_state(v, target, self.get_fold_sequence(self.sequence, target), colormap, filename)

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

    def get_design_score(self, sequences, energies=False):
        """
        calculates design score using scoring function
        """
        # in a small number of cases, get_design function causes error
        # if this happens, assume 0 score
        #try:
        if energies:
            return energies[0] - energies[1]
        if not self.scoring_func:
            return 0
        score = self.scoring_func(sequences)
        return score
        #except:
        #    return 0

    def get_sequence_info(self, sequence):
        """
        get sequence information - native fold, bp distance, pairmap, design
        for a particular sequence
        """
        native = [] 
        energies = []
        fold_sequences = []
        bpps = []
        if self.mode == "nupack":
            result = [None] * self.n_targets
            p = multiprocessing.Pool(self.n_targets)
        for i in range(self.n_targets):
            fold_sequence = self.get_fold_sequence(sequence, self.targets[i])
            fold_sequences.append(fold_sequence)
            if self.mode == "vienna":
                if self.targets[i]['type'] == "aptamer":
                    fold_list = fold_utils.vienna_fold(fold_sequences[i], self.targets[i]['fold_constraint'], bpp=True)
                else:
                    fold_list = fold_utils.vienna_fold(fold_sequences[i], bpp=True)
                native.append(fold_list[0])
                bpps.append(fold_list[2])
                if self.aptamer:
                    energies.append(fold_list[1])
            if self.mode == "nupack":
                if self.targets[i]['type'] == "oligos" and isinstance(self.targets[i]['inputs'], dict):
                    result[i] = p.apply_async(fold_utils.nupack_fold, args=(fold_sequence, [x*self.oligo_conc for x in self.targets[i]['inputs'].values()], True))
                else:
                    result[i] = p.apply_async(fold_utils.nupack_fold, args=(fold_sequence, self.target_oligo_conc*self.oligo_conc, True))
        if self.mode == "nupack":
            p.close()
            p.join()
            result = [x.get() for x in result]
            native = [[x[0], x[2]] for x in result]
            energies = [x[1] for x in result]
            bpps = [x[3] for x in result]
        if self.aptamer:
            design_score = self.get_design_score(bpps)
            bp_distance = self.score_secstructs(native, energies, sequence)
        else:
            design_score = max(self.get_design_score(bpps),0)
            bp_distance = self.score_secstructs(native, sequence=sequence)
        return [sequence, native, bp_distance, fold_sequences, design_score]

    def reset_sequence(self):
        """
        reset sequence to the start sequence (for rerunning optimization)
        """
        self.sequence_graph.reset_sequence(self.beginseq)
        self.update_current(self.beginseq)
        self.update_best()
        self.all_solutions = []
        if self.print_:
            print "reset %s" % self.sequence

    def update_current(self, sequence, native=None, bp_distance=None, design_score=None, energies=None):
        """
        updates current sequence and related information
        """
        self.sequence = sequence
        if not native:
            [sequence, native, bp_distance, fold_sequences, design_score] = self.get_sequence_info(sequence)
        self.native = native
        self.bp_distance = bp_distance
        self.design_score = design_score

    def update_best(self):
        """
        updates best to current sequence
        """
        self.best_sequence = self.sequence
        self.best_native = self.native
        self.best_bp_distance = self.bp_distance
        self.best_design_score = self.design_score

    def score_secstructs(self, secstruct, energies = False, sequence = False):
        """
        calculates sum of bp distances for with and without oligo

        returns:
        sum of bp distances with and without the oligo 
        """
        # test for secondary structure matches
        distance = 0.0
        strands_interacting = 0.0
        n_strands = 0.0
        for i in range(self.n_targets):
            if "threshold" in self.targets[i]:
                distance += self.bp_distance_func(secstruct[i], self.targets[i]['secstruct'], self.targets[i]['constrained'], self.targets[i]['threshold'])
            else:
                distance += self.bp_distance_func(secstruct[i], self.targets[i]['secstruct'], self.targets[i]['constrained'])
            if self.mode == "nupack":
                strands = secstruct[i][0].split("&")
                for strand in secstruct[i][1]:
                    if strand != len(strands):
                        n_strands += 1
                        if "(" in strands[strand-1] or ")" in strands[strand-1]:
                            strands_interacting += 1

        # add bonus for interaction of strands
        if self.strandbonus:
            if strands_interacting == 0:
                strands_interacting += 1
            distance /= strands_interacting/n_strands
            if self.print_:
                print "bonus: %d" % distance

        # test energies
        if energies:
            if energies[1] - 0.6 * math.log(self.aptamer/3.0) > energies[0]:
                distance += 4

        # test sequence
        if sequence:
            distance += 1 * sequence.count("GGGG")
            distance += 1 * sequence.count("CCCC")
            distance += 1 * sequence.count("AAAAA")

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
        return self.score_secstructs(self.best_native) == 0 and self.oligo_conc == 1

    def optimize_sequence(self, n_iterations, n_cool = 50, greedy = None, print_ = None, start_oligo_conc=1e7, continue_opt=False):
        """
        monte-carlo optimization of the sequence

        args:
        n_interations is the total number of iterations
        n_cool is the number of times to cool the system
        """
        bases = "GAUC"
        pairs = ["GC", "CG", "AU", "UA"]
    
        if greedy != None:
            self.greedy = greedy
        if print_ != None:
            self.print_ = print_

        #print self.targets
        #self.optimize_start_sequence()

        T = 5.0
        self.oligo_conc = start_oligo_conc

        def p_dist(dist, new_dist):
            """probability function"""
            return math.exp(-(new_dist-dist)/T)

        def p_greedy(dist, new_dist):
            if new_dist <= dist:
                return 1
            else:
                return 0

        if self.greedy:
            p_func = p_greedy
        else:
            p_func = p_dist

        niter = None

        # loop as long as bp distance too large or design score too small
        for i in range(n_iterations):
            #random.shuffle(index_array)
            
            # pick random nucleotide in sequence
            mut_sequence = self.sequence_graph.mutate()
            [mut_sequence, native, bp_distance, fold_sequences, design_score] = self.get_sequence_info(mut_sequence)

            # if distance or score is better for mutant, update the current sequence
            p = p_func(self.bp_distance, bp_distance)
            if(random.random() <= p):
                if self.bp_distance == bp_distance and p_func(self.design_score, design_score):
                    continue
                self.update_current(mut_sequence, native, bp_distance, design_score)
                if self.print_:
                    print self.sequence, self.bp_distance, self.design_score
                    print "conc: %s" % self.oligo_conc
                    for j in range(self.n_targets):
                        print self.native[j]
                    print ""
                        #print self.get_fold_sequence(self.sequence, self.targets[j])
            
                # if distance or score is better for mutant than best, update the best sequence    
                if(bp_distance < self.best_bp_distance or
                   (bp_distance == self.best_bp_distance and design_score > self.best_design_score)):
                    self.update_best()

            if self.best_bp_distance == 0 and self.oligo_conc == 1.0:
                print "-> Reached solution in %d iterations." % i
                if not continue_opt:
                    return i

            # decrease temperature
            #if i % (n_iterations/n_cool) == 0:
            wait = 0
            interval = n_iterations/(2*n_cool)
            if i % interval == 0 and i >= interval*wait and i < interval*(n_cool+wait):
                T -= 0.1
                if T < 1:
                    T = 1
            
            # update oligo_conc
            while self.best_bp_distance == 0 and self.oligo_conc != 1.0:
                if self.oligo_conc/10 <= 1.0:
                    self.oligo_conc = 1.0
                else:
                    self.oligo_conc /= 10
                self.update_current(self.sequence)
                self.update_best()
        
        return niter

