import design_utils, fold_utils
import sequence_graph
import settings
import math, random
import unittest
import sys
import multiprocessing

class SwitchDesigner(object):

    def __init__(self, id, beginseq, constraints, targets, **kwargs):
        # sequence information
        self.id = id
        self.beginseq = beginseq
        self.sequence = beginseq
        self.n = len(self.sequence)
        self.constraints = constraints

        # read keyword arguments
        self.mode = kwargs.get('mode', 'nupack')
        add_rcs = kwargs.get('add_rcs', False)
        self.print_ = kwargs.get('print_', False)
        self.inputs = kwargs.get('inputs', {})
        self.substrings = kwargs.get('substrings', [])

        # automatic folding mode detection
        aptamer = False
        if any([target['type'] == 'aptamer' for target in targets]):
            aptamer = True
        multi_oligos = False
        if any([len(target['inputs']) > 1 for target in targets if 'inputs' in target]):
            multi_oligos = True
        if aptamer and multi_oligos:
            print 'Unable to handle aptamers and multiple RNA inputs simultaneously'
            sys.exit()
        elif aptamer:
            if self.mode != 'vienna':
                print 'Switching to Vienna to handle aptamer'
                self.mode = 'vienna'
        elif multi_oligos:
            if self.mode != 'nupack':
                print 'Switching to NUPACK to handle multiple RNA inputs'
                self.mode = 'nupack'

        # process input data
        self.targets = self.parse_targets(targets)
        self.n_targets = len(self.targets)
        self.sequence_graph = sequence_graph.SequenceGraph(self.inputs, targets, constraints, beginseq, add_rcs, False, autocomplement=True)
        self.target_oligo_conc = 1e-7

        # initialize default optimization parameters
        self.scoring_func = design_utils.get_bpp_scoring_func(targets, self.mode == 'nupack')
        self.greedy = False
        self.oligo_conc = 1.0

        # print design info
        if self.print_:
            print self.constraints
            for i, target in enumerate(targets):
                print '-> state %d' % i
                print self.get_fold_sequence(self.sequence, target)
                print target['secstruct']
                print target['constrained']
                if 'fold_constraint' in target:
                    print target['fold_constraint']

        self.update_current(self.sequence)
        
        # maintain sequences
        self.update_best()

    def parse_targets(self, targets):
        """
        generate full secondary structure and constraint strings
        """
        for target in targets:
            secstruct = ''
            constrained = ''
            if target['type'] == 'oligos' and len(target['inputs']) > 0 and '&' not in target['secstruct']:
                for input in target['inputs']:
                    if self.inputs[input]['type'] == 'RNA':
                        secstruct += '.'*len(self.inputs[input]['sequence']) + '&'
                        constrained += 'o'*len(self.inputs[input]['sequence']) + 'x'
                target['secstruct'] = secstruct + target['secstruct']
                target['constrained'] = constrained + target['constrained']
            #elif target['type'] == 'aptamer':
            #    self.aptamer = float(target['concentration'])
            #    fold_constraint = list(target['secstruct'])
            #    for i, fold in enumerate(fold_constraint):
            #        if i in target['site'] and fold == '.':
            #            fold_constraint[i] = 'x'
            #        elif i not in target['site']:
            #            fold_constraint[i] = '.'
            #    target['fold_constraint'] = ''.join(fold_constraint)
            if '&' in target['secstruct']:
                breaks = [i for i, char in enumerate(target['secstruct']) if char == '&']
                constrained = list(target['constrained'])
                for i in breaks:
                    constrained[i] = 'x'
                target['constrained'] = ''.join(constrained)
        return targets

    def get_fold_sequence(self, sequence, target):
        """ append oligo sequences separated by & for type oligo """
        if target['type'] == 'oligos':
            return '&'.join([self.inputs[x]['sequence'] for x in sorted(target['inputs']) if self.inputs[x]['type'] == 'RNA' ] + [sequence])
        else:
            return sequence 

    def get_solution(self):
        """
        return current best as solution
        """
        for i, target in enumerate(self.targets):
            fold_sequence = self.get_fold_sequence(self.best_sequence, target)
            print fold_sequence
            if self.mode == 'vienna':
                if target['type'] == 'aptamer':
                    ligand = self.inputs[target['inputs'].keys()[0]]
                    print fold_utils.vienna_fold(fold_sequence, ligand['fold_constraint'], bpp=False)[0]
                else:
                    print fold_utils.vienna_fold(fold_sequence)[0]
            else:
                concentrations = [target['inputs'][input]*self.oligo_conc for input in sorted(target['inputs'])]
                print fold_utils.nupack_fold(fold_sequence, concentrations, False)[0]
        return [self.best_sequence, self.best_bp_distance, self.best_design_score]

    def get_bpp_score(self, bpps):
        """
        calculates design score using scoring function
        """
        score = self.scoring_func(bpps)
        return score

    def get_sequence_info(self, sequence):
        """
        get sequence information - native fold, bp distance, pairmap, design
        for a particular sequence
        """
        native = [] 
        energies = []
        fold_sequences = []
        bpps = []
        if self.mode == 'nupack':
            result = [None] * self.n_targets
            p = multiprocessing.Pool(self.n_targets)
        for i, target in enumerate(self.targets):
            fold_sequence = self.get_fold_sequence(sequence, target)
            fold_sequences.append(fold_sequence)
            if self.mode == 'vienna':
                if target['type'] == 'aptamer':
                    ligand = self.inputs[target['inputs'].keys()[0]]
                    fold_list = fold_utils.vienna_fold(fold_sequences[i], ligand['fold_constraint'], bpp=True)
                else:
                    fold_list = fold_utils.vienna_fold(fold_sequences[i], bpp=True)
                native.append(fold_list[0])
                energies.append(fold_list[1])
                bpps.append(fold_list[2])
            if self.mode == 'nupack':
                #if self.targets[i]['type'] == 'oligos' and isinstance(self.targets[i]['inputs'], dict):
                concentrations = [target['inputs'][input]*self.oligo_conc for input in sorted(target['inputs'])]
                result[i] = p.apply_async(fold_utils.nupack_fold, args=(fold_sequence, concentrations, True))
                #else:
                #    result[i] = p.apply_async(fold_utils.nupack_fold, args=(fold_sequence, self.target_oligo_conc*self.oligo_conc, True))
        if self.mode == 'nupack':
            p.close()
            p.join()
            result = [x.get() for x in result]
            native = [[x[0], x[2]] for x in result]
            energies = [x[1] for x in result]
            bpps = [x[3] for x in result]
        design_score = max(self.get_bpp_score(bpps),0)
        bp_distance = self.score_secstructs(sequence, native, energies)
        return [sequence, native, energies, bp_distance, fold_sequences, design_score]

    def reset_sequence(self):
        """
        reset sequence to the start sequence (for rerunning optimization)
        """
        self.sequence_graph.reset_sequence(self.beginseq)
        self.update_current(self.beginseq)
        self.update_best()
        if self.print_:
            print 'reset %s' % self.sequence

    def update_current(self, sequence, native=None, energies=None, bp_distance=None, design_score=None):
        """
        updates current sequence and related information
        """
        self.sequence = sequence
        if not native:
            [sequence, native, energies, bp_distance, fold_sequences, design_score] = self.get_sequence_info(sequence)
        self.native = native
        self.energies = energies
        self.bp_distance = bp_distance
        self.design_score = design_score

    def update_best(self):
        """
        updates best to current sequence
        """
        self.best_sequence = self.sequence
        self.best_native = self.native
        self.best_energies = self.energies
        self.best_bp_distance = self.bp_distance
        self.best_design_score = self.design_score

    def score_secstructs(self, sequence, secstruct, energies):
        """
        calculates sum of bp distances for with and without oligo

        returns:
        sum of bp distances with and without the oligo 
        """
       # test for secondary structure matches
        energy_compare = {}
        distance = 0.0
        strands_interacting = 0.0
        n_strands = 0.0
        for i, target in enumerate(self.targets):
            if 'threshold' in target:
                distance += design_utils.bp_distance(secstruct[i], target['secstruct'], target['constrained'], target['threshold'])
            else:
                distance += design_utils.bp_distance(secstruct[i], target['secstruct'], target['constrained'])
            if target['type'] == 'aptamer':
                energy_compare[target['type']] = energies[i]
                ligand = target['inputs'].keys()[0]
                energy_compare['ligand'] = self.inputs[ligand]['kD'], target['inputs'][ligand]
            elif target['type'] == 'single':
                energy_compare[target['type']] = energies[i]

        # test energies
        if 'aptamer' in energy_compare:
            if energy_compare['aptamer'] - 0.6 * math.log(energy_compare['ligand'][1]/3.0) > energy_compare['single']:
                distance += 4

        for substr in self.substrings:
            if substr in sequence:
                distance += 1

        return distance

    def check_current_secstructs(self):
        return self.score_secstructs(self.best_sequence, self.best_native, self.best_energies) == 0 and self.oligo_conc == 1

    def optimize_sequence(self, n_iterations, n_cool = 50, greedy = None, print_ = None, start_oligo_conc=1e7, continue_opt=False):
        """
        monte-carlo optimization of the sequence

        args:
        n_interations is the total number of iterations
        n_cool is the number of times to cool the system
        """
        bases = 'GAUC'
        pairs = ['GC', 'CG', 'AU', 'UA']
    
        if greedy != None:
            self.greedy = greedy
        if print_ != None:
            self.print_ = print_

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
            [mut_sequence, native, energies, bp_distance, fold_sequences, design_score] = self.get_sequence_info(mut_sequence)

            # if distance or score is better for mutant, update the current sequence
            p = p_func(self.bp_distance, bp_distance)
            if(random.random() <= p):
                if self.bp_distance == bp_distance and random.random() > p_func(design_score, self.design_score):
                    continue
                self.update_current(mut_sequence, native, energies, bp_distance, design_score)
                if self.print_:
                    print self.sequence, self.bp_distance, self.design_score
                    print 'conc: %s' % self.oligo_conc
                    for j in range(self.n_targets):
                        print self.native[j]
                    print ''
                        #print self.get_fold_sequence(self.sequence, self.targets[j])
            
                # if distance or score is better for mutant than best, update the best sequence    
                if(bp_distance < self.best_bp_distance or
                   (bp_distance == self.best_bp_distance and design_score > self.best_design_score)):
                    self.update_best()

            if self.best_bp_distance == 0 and self.oligo_conc == 1.0:
                print '-> Reached solution in %d iterations.' % i
                if not continue_opt:
                    return i

            # decrease temperature
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

