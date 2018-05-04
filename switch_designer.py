import design_utils, fold_utils
import sequence_graph
import settings
import math, random
import unittest
import sys
import multiprocessing

class Design(object):
    
    def __init__(self, begin_seq, seq_locks, targets, inputs, **kwargs):
        self.substrings = kwargs.get('substrings', [])
        
        self.begin_seq = begin_seq
        self.seq_locks = seq_locks

        self.inputs = inputs
        self.targets = self.parse_targets(targets)
        self.n_targets = len(self.targets)
        self.get_default_mode()

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

    def get_fold_sequences(self, sequence):
        """ append oligo sequences separated by & for type oligo """
        fold_sequences = []
        for target in self.targets:
            if target['type'] == 'oligos':
                fold_sequences.append('&'.join([self.inputs[x]['sequence'] for x in sorted(target['inputs']) if self.inputs[x]['type'] == 'RNA' ] + [sequence]))
            else:
                fold_sequences.append(sequence)
        return fold_sequences
    
    def get_default_mode(self):
        aptamer = False
        if any([target['type'] == 'aptamer' for target in self.targets]):
            aptamer = True
        multi_oligos = False
        if any([len(target['inputs']) > 1 for target in self.targets if 'inputs' in target]):
            multi_oligos = True
        if multi_oligos:
            print 'Switching to NUPACK to handle multiple RNA inputs'
            self.default_mode = 'nupack'
        else:
            self.default_mode = False

        return self.default_mode

    def print_info(self):        
        print self.seq_locks
        fold_sequences = self.get_fold_sequences(self.begin_seq)
        for i, target in enumerate(self.targets):
            print '-> state %d' % i
            print fold_sequences[i]
            print target['secstruct']
            print target['constrained']
            if 'fold_constraint' in target:
                print target['fold_constraint']

class DesignSequence(object):

    def __init__(self, design, sequence, mode = 'nupack', oligo_conc = 1):
        self.design = design
        self.n_targets = len(design.targets)
        self.mode = design.default_mode if design.default_mode else mode
        self.scoring_func = design_utils.get_bpp_scoring_func(design.targets, self.mode == 'nupack')

        self.update_sequence(sequence, oligo_conc)

    def get_design_score(self):
        """
        calculates design score using scoring function
        """
        return self.scoring_func(self)
    
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
        self.mispaired_positions = set()
        for i, target in enumerate(self.design.targets):
            if 'threshold' in target:
                bp_result = design_utils.bp_distance(secstruct[i], target['secstruct'], target['constrained'], target['threshold'])
                distance += bp_result[0]
                offset = len(target['secstruct']) - len(sequence)
                self.mispaired_positions.update(set([j - offset for j in bp_result[1] if j >= offset]))
            else:
                bp_result = design_utils.bp_distance(secstruct[i], target['secstruct'], target['constrained'])
                distance += bp_result[0]
                offset = len(target['secstruct']) - len(sequence)
                self.mispaired_positions.update(set([j - offset for j in bp_result[1] if j >= offset]))
            if target['type'] == 'aptamer':
                energy_compare[target['type']] = energies[i]
                ligand = target['inputs'].keys()[0]
                energy_compare['ligand'] = self.design.inputs[ligand]['kD'], target['inputs'][ligand]
            elif target['type'] == 'single':
                energy_compare[target['type']] = energies[i]

        # test energies
        if 'aptamer' in energy_compare:
            if energy_compare['aptamer'] - 0.6 * math.log(energy_compare['ligand'][1]/energy_compare['ligand'][0]) > energy_compare['single']:
                distance += 4

        for substr in self.design.substrings:
            if substr in sequence:
                distance += 1

        return distance

    def update_sequence(self, sequence, oligo_conc=1):
        self.sequence = sequence
        self.native = [] 
        self.energies = []
        self.fold_sequences = self.design.get_fold_sequences(sequence)
        self.bpps = []
        self.oligo_conc = oligo_conc
        if self.mode == 'nupack':
            result = [None] * self.n_targets
            p = multiprocessing.Pool(self.n_targets)
        for i, target in enumerate(self.design.targets):
            if self.mode == 'vienna':
                if target['type'] == 'aptamer':
                    ligand = self.design.inputs[target['inputs'].keys()[0]]
                    fold_list = fold_utils.vienna_fold(self.fold_sequences[i], ligand['fold_constraint'], bpp=True)
                else:
                    fold_list = fold_utils.vienna_fold(self.fold_sequences[i], bpp=True)
                self.native.append(fold_list[0])
                self.energies.append(fold_list[1])
                self.bpps.append(fold_list[2])
            if self.mode == 'nupack':
                if 'inputs' in target:
                    concentrations = [target['inputs'][input]*oligo_conc for input in sorted(target['inputs'])]
                else:
                    concentrations = 1
                if target['type'] == 'aptamer':
                    ligand = self.design.inputs[target['inputs'].keys()[0]]
                    result[i] = p.apply_async(fold_utils.nupack_fold,
                                              args=(self.fold_sequences[i],
                                                    ligand['fold_constraint'],
                                                    concentrations, True))
                else:
                    result[i] = p.apply_async(fold_utils.nupack_fold,
                                              args=(self.fold_sequences[i], False,
                                                    concentrations, True))
        if self.mode == 'nupack':
            p.close()
            p.join()
            result = [x.get() for x in result]
            self.native = [[x[0], x[2]] for x in result]
            self.energies = [x[1] for x in result]
            self.bpps = [x[3] for x in result]
        self.bp_distance = self.score_secstructs(sequence, self.native, self.energies)
        self.design_score = max(self.get_design_score(),0)
        return
    
    def is_solution(self):
        return self.score_secstructs(self.sequence, self.native, self.energies) == 0 and self.oligo_conc == 1

    def print_(self):
        print self.sequence
        print 'bp distance: %d' % self.bp_distance
        print 'design score: %f' % self.design_score
        print 'conc: %s' % self.oligo_conc
        for j in range(self.n_targets):
            print self.native[j]
        print ''

class SwitchDesigner(object):

    def __init__(self, id, design, **kwargs):
        self.id = id
        self.design = design

        # read keyword arguments
        self.mode = kwargs.get('mode', 'nupack')
        add_rcs = kwargs.get('add_rcs', False)
        self.print_ = kwargs.get('print_', False)
        self.inputs = kwargs.get('inputs', {})

        # process input data
        self.sequence_graph = sequence_graph.SequenceGraph(self.design, add_rcs=add_rcs)
        self.target_oligo_conc = 1e-7

        # initialize default optimization parameters
        self.greedy = False
        self.oligo_conc = 1.0

        if self.print_:
            self.design.print_info()

        # set designs to start
        self.reset_sequence()


    def get_solution(self):
        """
        return current best as solution
        """
        return self.current_design

    def reset_sequence(self):
        """
        reset sequence to the start sequence (for rerunning optimization)
        """
        sequence = self.design.begin_seq
        self.sequence_graph.reset_sequence(sequence)
        self.current_design = DesignSequence(self.design, sequence, self.mode, self.oligo_conc)
        self.update_best()
        if self.print_:
            print 'reset %s' % sequence

    def update_current(self, design):
        """
        updates current sequence and related information
        """
        self.current_design = design

    def update_best(self):
        """
        updates best to current sequence
        """
        self.best_design = self.current_design

    def check_current(self):
        return self.best_design.is_solution()

    def optimize_sequence(self, n_iterations, n_cool = 50, greedy = None, print_ = None, start_oligo_conc=1e7, continue_=False):
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
            mut_sequence = self.sequence_graph.mutate(self.current_design.mispaired_positions)
            mut_design = DesignSequence(self.design, mut_sequence, self.mode, self.oligo_conc)

            # if distance or score is better for mutant, update the current sequence
            p = p_func(self.current_design.bp_distance, mut_design.bp_distance)
            if(random.random() <= p):
                if self.current_design.bp_distance == mut_design.bp_distance and random.random() > p_func(mut_design.design_score, self.current_design.design_score):
                    continue
                self.update_current(mut_design)
                if self.print_:
                    print i
                    self.current_design.print_()
                        #print self.get_fold_sequence(self.sequence, self.targets[j])
            
                # if distance or score is better for mutant than best, update the best sequence    
                if(self.current_design.bp_distance < self.best_design.bp_distance or
                   (self.current_design.bp_distance == self.best_design.bp_distance and 
                    self.current_design.design_score > self.best_design.design_score)):
                    self.update_best()

            if self.best_design.bp_distance == 0 and self.oligo_conc == 1.0:
                print '-> Reached solution in %d iterations.' % i
                self.best_design.print_()
                if not continue_:
                    return i

            # decrease temperature
            wait = 0
            interval = n_iterations/(2*n_cool)
            if i % interval == 0 and i >= interval*wait and i < interval*(n_cool+wait):
                T -= 0.1
                if T < 1:
                    T = 1
            
            # update oligo_conc
            while self.best_design.bp_distance == 0 and self.oligo_conc != 1.0:
                if self.oligo_conc/10 <= 1.0:
                    self.oligo_conc = 1.0
                else:
                    self.oligo_conc /= 10
                self.current_design.update_sequence(self.sequence, self.oligo_conc)
                self.update_best()
        
        return niter

