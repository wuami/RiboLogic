import eterna_utils
import design_utils
import varna, draw_utils
import settings
import inv_utils
import random
import math
import itertools
import ensemble_design
import unittest
import sys
import networkx as nx
import multiprocessing

class SwitchDesigner(object):

    def __init__(self, id, type, beginseq, constraints, targets, design_linker, scoring = "bpp", inputs = None, mode = "ghost", oligorc = False):
        # sequence information
        self.id = id
        self.type = type
        self.beginseq = beginseq
        self.sequence = beginseq
        self.n = len(self.sequence)
        self.constraints = constraints
        self.index_array = self.get_unconstrained_indices()

        # scoring
        if scoring == "bpp":
            self.scoring_func = design_utils.get_bpp_scoring_func(targets)
        elif scoring == "ensemble":
            self.scoring_func = design_utils.get_ensemble_scoring_func()
        elif scoring == "landing":
            self.scoring_func = design_utils.get_strategy_scoring_func("eli_landing_lane")
        elif scoring == "null":
            self.scoring_func = None
        else:
            raise ValueError("invalid scoring function: %s" % scoring)

        #if type == "multi_input_oligo":
        self.bp_distance_func = design_utils.bp_distance_with_unpaired
        self.greedy = False
        self.oligo_conc = 1e-7
        #else:
        #    self.bp_distance_func = design_utils.bp_distance_with_constraint
        self.mode = mode
        if self.mode == "hairpin":
            self.hp_mismatch = False

        # target information
        self.targets = targets
        self.n_targets = len(self.targets)
        self.inputs = inputs
        if "pos" in self.inputs:
            self.input_pos = self.inputs['pos']
            self.input_pos.insert(0,0)
            del self.inputs['pos']
        else:
            self.input_pos = [0]*(len(self.inputs)+1)
        n_cores = min(16, self.n_targets)
        
        # update dependency graph
        self.dep_graph = self.get_dependency_graph()
        seq_array = ensemble_design.get_sequence_array(self.sequence)
        for i in range(self.n):
            if self.constraints[i] == "x":
                self.update_neighbors(i, seq_array, [])
        self.sequence = ensemble_design.get_sequence_string(self.sequence)

        if mode == "nupack":
            self.linker_length = 0
            self.linker = ""
            self.design_linker = ""
        else:
            self.linker_length = 5
            self.linker = ensemble_design.get_sequence_string(["C" if i % 5 == 2 else "A" for i in range(self.linker_length)])
            self.free_linkers = True
            self.design_linker = design_linker
            self.oligotail = ""
        if type == "multi_input" or type == "multi_input_oligo":
            self.create_target_secstructs()
        self.cotrans = False
        
        # print puzzle info
        self.prints = False
        for i in targets:
            print self.get_fold_sequence(self.sequence, i)
            print i['secstruct']
            print i['constrained']

        self.update_sequence(*self.get_sequence_info(self.sequence))
        
        # maintain sequences
        self.update_best()
        self.all_solutions = []

        #if type == "multi_input_oligo":
        if oligorc:
            self.set_oligo_rcs()
            self.mutate_func = self.mutate_or_shift
        else:
            self.mutate_func = self.mutate_sequence
        
    def get_dependency_graph(self):
        """ get dependency graph based on target secondary structures """
        graph = nx.Graph()
        graph.add_nodes_from(range(self.n))
        for target in self.targets:
            pairmap = eterna_utils.get_pairmap_from_secstruct(target['secstruct'].split('&')[-1])
            for i, j in enumerate(pairmap):
                if j > i:
                    graph.add_edge(i,j)
        if not nx.is_bipartite(graph):
            raise ValueError("dependency graph is not bipartite")
        return graph

    def set_oligo_rcs(self):
        self.oligo_rc = []
        self.oligo_pos = []
        self.oligo_len = []
        self.oligo_len_sum = 0
        inputs = self.inputs.values()
        self.set_oligo_rc(inputs[0], [0,len(inputs[0])])
        self.set_oligo_rc(inputs[1], [self.n-len(inputs[1]), self.n])

    def set_oligo_rc(self, seq, range):
        """ set oligo rc parameters for shifting complement sequece"""
        #lo, hi = self.n-49, self.n-27
        rc = design_utils.rc(seq)
        self.beginseq = self.beginseq[0:range[0]] + rc + self.beginseq[range[1]:]
        self.oligo_rc.append(rc)
        self.oligo_pos.append(range)
        self.oligo_len.append([0,range[1]-range[0]-1])
        self.oligo_len_sum += len(seq)

    def create_target_secstructs(self):
        """ add oligo secstructs to target secstruct """
        for target in self.targets:
            if target['type'] == 'oligos' and target['secstruct'] == self.n:
                # initialize sequences
                secstruct = ""
                constrained = ""
                if self.mode == "ghost":
                    fold_constraint = ""
                # get positions of inputs
                for i,input in enumerate(sorted(self.inputs)):
                    # add portion of sequence between inputs
                    secstruct += target['secstruct'][self.input_pos[i]:self.input_pos[i+1]]
                    constrained += target['constrained'][self.input_pos[i]:self.input_pos[i+1]]
                    if self.mode == "ghost":
                        fold_constraint += '.'*(self.input_pos[i+1]-self.input_pos[i])
                    if self.input_pos[i+1] - self.input_pos[i] != 0:
                        secstruct += '.'*len(self.design_linker)
                        if self.free_linkers:
                            constrained += 'o'*len(self.design_linker)
                        else:
                            constrained += 'u'*len(self.design_linker)
                        if self.mode == "ghost":
                            fold_constraint += '.'*len(self.design_linker)
                    # add input sequence
                    seq = self.inputs[input]
                    if input in target['inputs']:
                        constrained += 'o'*len(seq)
                        if self.mode == "hairpin":
                            hairpin_len = (len(seq)-4)/2
                            if self.hp_mismatch:
                                mismatch = hairpin_len/2
                                secstruct += '('*(mismatch-1)+'.'+'('*(hairpin_len-mismatch) + '.'*(len(seq)-2*hairpin_len) + ')'*(hairpin_len-mismatch) + '.' + ')'*(mismatch-1)
                            else:
                                secstruct += '('*hairpin_len + '.'*(len(seq)-2*hairpin_len) + ')'*hairpin_len
                        elif self.mode == "ghost":
                            secstruct += '.'*len(seq)
                            fold_constraint += 'x'*len(seq)
                        elif self.mode == "nupack":
                            secstruct += '.'*len(seq) + '&'
                            constrained += '&'
                    elif self.mode != "nupack":
                        secstruct += '.'*len(seq)
                        constrained += 'o'*len(seq)
                        if self.mode == "ghost":
                            fold_constraint += '.'*len(seq)
                    # add linker sequence
                    if self.mode != "nupack" and i != len(self.inputs)-1:
                        secstruct += '.'*self.linker_length
                        if self.free_linkers:
                            constrained += 'o'*self.linker_length
                        else:
                            constrained += 'u'*self.linker_length
                        if self.mode == "ghost":
                            fold_constraint += '.'*self.linker_length
                if self.mode != "nupack" and self.input_pos[-1] != len(target['secstruct']):
                    secstruct += '.'*len(self.design_linker)
                    if self.free_linkers:
                        constrained += 'o'*len(self.design_linker)
                    else:
                        constrained += 'u'*len(self.design_linker)
                    if self.mode == "ghost":
                        fold_constraint += '.'*len(self.design_linker)
                # add last portion of sequence
                target['secstruct'] = secstruct + target['secstruct'][self.input_pos[-1]:]
                target['constrained'] = constrained + target['constrained'][self.input_pos[-1]:]
                if self.mode == "ghost":
                    target['fold_constraint'] = fold_constraint + '.'*(self.n-self.input_pos[-1])

    def get_fold_sequence(self, sequence, objective):
        """ append oligo sequences separated by & for type oligo """
        if objective['type'] == 'oligo':
            return '&'.join([sequence, objective['oligo_sequence']])
        elif objective['type'] == 'oligos':
            if self.mode == "nupack":
                return '&'.join([self.inputs[x] for x in objective['inputs']] + [sequence])
            # get positions of inputs
            fold_seq = ""
            for i,input in enumerate(sorted(self.inputs)):
                # add portion of sequence between two inputs
                fold_seq += sequence[self.input_pos[i]:self.input_pos[i+1]]
                if self.input_pos[i+1]-self.input_pos[i] != 0:
                    fold_seq += self.design_linker
                # add input sequence
                if input in objective['inputs']:
                    if self.mode == "ghost":
                        fold_seq += design_utils.rc(self.inputs[input])
                    elif i == len(self.inputs)-1:
                        fold_seq += self.get_hairpin(design_utils.rc(self.inputs[input]), False, self.hp_mismatch)
                    else:
                        fold_seq += self.get_hairpin(design_utils.rc(self.inputs[input]), True, self.hp_mismatch)
                else:
                    fold_seq += design_utils.rc(self.inputs[input])
                # add linker sequence
                if i != len(self.inputs)-1:
                    fold_seq += self.linker
            if self.input_pos[-1] != len(sequence):
                fold_seq += self.design_linker
            return fold_seq + sequence[self.input_pos[-1]:]
        else:
            return sequence 

    def get_unconstrained_indices(self):
        """
        get indices for positions that can change)
        """
        open = []
        restricted = []
        for ii in range(0,self.n):
            if(self.constraints[ii] == "o"):
                open.append(ii)
            elif self.constraints[ii] == "r":
                restricted.append(ii)
        return [open, restricted]

    def get_solution(self):
        """
        return current best as solution
        """
        for i in range(self.n_targets):
            fold_sequence = self.get_fold_sequence(self.best_sequence, self.targets[i])
            print fold_sequence
            if self.mode == "ghost":
                print inv_utils.fold(fold_sequence, self.cotrans, self.targets[i]['fold_constraint'])[0]
            elif self.mode == "hairpin":
                print inv_utils.fold(fold_sequence, self.cotrans)[0]
            else:
                print inv_utils.nupack_fold(fold_sequence, self.target_oligo_conc)[0]
        return [self.best_sequence, self.best_bp_distance, self.best_design_score]

    def draw_solution(self, name):
        """ draw each state for current solution """
        # initiate varna RNA visualizer
        v = varna.Varna()
        
        # get puzzle object and generate colormaps for each objective
        n = len(self.inputs)
        colormap = draw_utils.get_colormaps(self.targets, self.inputs, self.input_pos, self.n, self.linker_length, self.design_linker, n)
        
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

    def get_design_score(self, sequence, secstruct):
        """
        calculates design score using scoring function
        """
        # in a small number of cases, get_design function causes error
        # if this happens, assume 0 score
        try:
            if not self.scoring_func:
                return 0
            designs = []
            for i in range(self.n_targets):
                if self.mode == "ghost":
                    designs.append(eterna_utils.get_design_from_sequence(sequence[i], secstruct[i], self.targets[i]['fold_constraint']))
                elif self.mode == "nupack":
                    designs.append(eterna_utils.get_design_from_sequence(sequence[i], secstruct[i], nupack=True))
                else:
                    designs.append(eterna_utils.get_design_from_sequence(sequence[i], secstruct[i]))
            score = self.scoring_func(designs)
            return score
        except:
            return 0

    def get_hairpin(self, seq, begin, mismatch):
        """ turn sequence into hairpin """
        n = (len(seq)-4)/2
        if begin:
            hairpin = seq[0:n] + seq[n:len(seq)-n] + design_utils.rc(seq[0:n])
            i = -n/2
        else:
            hairpin = design_utils.rc(seq[-n:]) + seq[n:len(seq)-n] + seq[-n:]
            i = n/2
        if mismatch:
            hairpin = ensemble_design.get_sequence_array(hairpin)
            hairpin[i] = design_utils.get_different_base(hairpin[i])
            return ensemble_design.get_sequence_string(hairpin)
        return hairpin

    def get_sequence_info(self, sequence):
        """
        get sequence information - native fold, bp distance, pairmap, design
        for a particular sequence
        """
        native = [""] * self.n_targets
        fold_sequences = []
        if self.mode == "nupack":
            p = multiprocessing.Pool(self.n_targets)
        for i in range(self.n_targets):
            fold_sequence = self.get_fold_sequence(sequence, self.targets[i])
            fold_sequences.append(fold_sequence)
            if self.mode == "ghost":
                fold = inv_utils.fold(fold_sequence, self.cotrans, self.targets[i]['fold_constraint'])[0]
                native.append(fold)
            elif self.mode == "hairpin":
                fold = inv_utils.fold(fold_sequence, self.cotrans)[0]
                native.append(fold)
            else:
                native[i] = p.apply_async(inv_utils.nupack_fold, args=(fold_sequence, self.oligo_conc))
        if self.mode == "nupack":
            p.close()
            p.join()
            native = [x.get() for x in native]
            native = [[x[0], x[2]] for x in native]
        bp_distance = self.score_secstructs(native)
        return [sequence, native, bp_distance, fold_sequences]

    def reset_sequence(self):
        """
        reset sequence to the start sequence (for rerunning optimization)
        """
        self.sequence = self.beginseq
        self.update_sequence(self.sequence)
        self.update_best()
        self.all_solutions = []
#        if self.type == "multi_input_oligo":
        if self.mode == "nupack":
            self.set_oligo_rcs()

    def update_sequence(self, sequence, native=None, bp_distance=None, score=None):
        """
        updates current sequence and related information
        """
        self.sequence = sequence
        if not native:
            [sequence, native, bp_distance, fold_sequences] = self.get_sequence_info(sequence)
            score = self.get_design_score(fold_sequences, native)
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
        distance = 0.0
        strands_bonus = 0.0
        n_strands = 0.0
        for i in range(self.n_targets):
            if "threshold" in self.targets[i]:
                distance += self.bp_distance_func(secstruct[i], self.targets[i]['secstruct'], self.targets[i]['constrained'], self.targets[i]['threshold'])
            else:
                distance += self.bp_distance_func(secstruct[i], self.targets[i]['secstruct'], self.targets[i]['constrained'])
            if self.mode == "nupack":
                strands = secstruct[i][0].split("&")
                for strand in secstruct[i][1]:
                    n_strands += 1
                    if strand != len(strands) and "(" in strands[strand]:
                        strands_bonus += 1
            if self.prints:
                print distance,
        if strands_bonus:
            distance /= strands_bonus/n_strands
            if self.prints:
                print "bonus: %d" % distance
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
        return self.score_secstructs(self.best_native) == 0 and self.oligo_conc == self.target_oligo_conc

    def mutate_or_shift(self, sequence):
        """
        either mutate sequence or shift the complement to the oligo in the design sequence
        """
        mut_array = ensemble_design.get_sequence_array(self.sequence)
        # mutate randomly wp 0.5, otherwise mutate oligo rc
        if (random.random() > 0.8): #float(self.oligo_len_sum)/sum([len(x) for x in self.index_array])):
            return self.mutate_sequence(sequence)
        else:
            # randomly choose an oligo
            roligo = random.getrandbits(1)

            # choose left or right
            maxed_out = False
            min_index = min(list(itertools.chain(*self.index_array)))
            max_index = max(list(itertools.chain(*self.index_array)))
            if self.oligo_pos[roligo][1] >= self.n or self.oligo_pos[roligo][1] >= max_index or \
               self.oligo_len[roligo][1] >= len(self.oligo_rc[roligo])-1 or self.constraints[self.oligo_pos[roligo][1]] == 'x':
                rindex = 0
                maxed_out = True
            if self.oligo_pos[roligo][0] <= 0 or self.oligo_pos[roligo][0] <= min_index or \
                 self.oligo_len[roligo][0] <= 0 or self.constraints[self.oligo_pos[roligo][0]-1] == 'x':
                rindex = 1
                maxed_out = maxed_out and True
            print rindex
            if 'rindex' not in locals():
                rindex = random.getrandbits(1)
            # wp 0.5 expand
            if (random.random() < 0.6 or self.oligo_len[roligo][1]-self.oligo_len[roligo][0] <= len(self.oligo_rc)-5) and \
                self.oligo_len[roligo][1]-self.oligo_len[roligo][0] != len(self.oligo_rc)-1 and not maxed_out:
                if rindex:
                    mut_array[self.oligo_pos[roligo][rindex]] = self.oligo_rc[roligo][self.oligo_len[roligo][rindex]]
                    self.oligo_pos[roligo][rindex] += 1
                    self.oligo_len[roligo][rindex] += 1
                else:
                    self.oligo_pos[roligo][rindex] -= 1
                    self.oligo_len[roligo][rindex] -= 1
                    #print self.oligo_pos[roligo][rindex]
                    #print self.oligo_len[roligo][rindex]
                    mut_array[self.oligo_pos[roligo][rindex]] = self.oligo_rc[roligo][self.oligo_len[roligo][rindex]]
            # otherwise shrink
            else:
                if rindex:
                    self.oligo_pos[roligo][rindex] -= 1
                    self.oligo_len[roligo][rindex] -= 1
                    mut_array[self.oligo_pos[roligo][rindex]] = ensemble_design.get_random_base()
                else:
                    #print self.oligo_pos[roligo]
                    #print self.oligo_len[roligo]
                    mut_array[self.oligo_pos[roligo][rindex]] = ensemble_design.get_random_base()
                    self.oligo_pos[roligo][rindex] += 1
                    self.oligo_len[roligo][rindex] += 1 
            self.oligo_len_sum = sum([x[1]-x[0] for x in self.oligo_len])
        return ensemble_design.get_sequence_string(mut_array)

    def mutate_sequence(self, sequence):
        """ mutate one random position """
        while True:
            mut_array = ensemble_design.get_sequence_array(self.sequence)
            if random.random() < 0.9 or len(self.index_array[1]) == 0 and len(self.index_array[0]) != 0:
                rindex = self.index_array[0][int(random.random() * len(self.index_array[0]))]
            else:
                rindex = self.index_array[1][int(random.random() * len(self.index_array[1]))]
            rbase = ensemble_design.get_random_base()
            mut_array[rindex] = rbase
            self.update_neighbors(rindex, mut_array, [])
            if design_utils.satisfies_constraints(mut_array, self.beginseq, self.constraints):
                break
        return ensemble_design.get_sequence_string(mut_array)

    def update_neighbors(self, node, mut_array, updated):
        complement = design_utils.rc(mut_array[node])
        for pos in nx.all_neighbors(self.dep_graph, node):
            mut_array[pos] = complement
            if pos not in updated:
                updated.append(pos)
                self.update_neighbors(pos, mut_array, updated)

    def optimize_sequence(self, n_iterations, n_cool = 50, greedy = None, cotrans = None, prints = None, target_oligo_conc=1e-7):
        """
        monte-carlo optimization of the sequence

        args:
        n_interations is the total number of iterations
        n_cool is the number of times to cool the system
        """
        bases = "GAUC"
        pairs = ["GC", "CG", "AU", "UA"]
    
        if len(self.index_array[0]) == 0 and len(self.index_array[1]) == 0:
            return

        if greedy != None:
            self.greedy = greedy
        if cotrans != None:
            self.cotrans = cotrans
        if prints != None:
            self.prints = prints

        #print self.targets
        #self.optimize_start_sequence()

        T = 5.0
        self.target_oligo_conc = target_oligo_conc
        self.oligo_conc = 1.0

        def p_dist(dist, new_dist):
            """probability function"""
            return math.exp(-(new_dist-dist)/T)

        def p_greedy(dist, new_dist):
            if new_dist <= dist:
                return 1
            else:
                return 0

        # loop as long as bp distance too large or design score too small
        for i in range(n_iterations):
            #random.shuffle(index_array)
            
            # pick random nucleotide in sequence
            mut_sequence = self.mutate_func(self.sequence)
            [mut_sequence, native, bp_distance, fold_sequences] = self.get_sequence_info(mut_sequence)

            if self.best_bp_distance != 0 and bp_distance == 0:
                print i
            
            # if distance or score is better for mutant, update the current sequence
            if self.greedy:
                p = p_greedy(self.bp_distance, bp_distance)
            else:
                p = p_dist(self.bp_distance, bp_distance)
            if(random.random() <= p):
                score = max(self.get_design_score(fold_sequences, native),0)
                if self.bp_distance == bp_distance and self.design_score != 0 and random.random() > p_dist(self.design_score, score):
                    continue
                self.update_sequence(mut_sequence, native, bp_distance, score)
                print self.oligo_pos
                if self.prints:
                    print self.sequence, self.bp_distance, self.design_score
                    print "conc: %s" % self.oligo_conc
                    for j in range(self.n_targets):
                        print self.native[j]
                        #print self.get_fold_sequence(self.sequence, self.targets[j])
            
                # if distance or score is better for mutant than best, update the best sequence    
                if(bp_distance < self.best_bp_distance or
                   (bp_distance == self.best_bp_distance and score > self.best_design_score)):
                    self.update_best()

            # decrease temperature
            #if i % (n_iterations/n_cool) == 0:
            wait = 0
            interval = n_iterations/(2*n_cool)
            if i % interval == 0 and i >= interval*wait and i < interval*(n_cool+wait):
                T -= 0.1
                if T < 1:
                    T = 1
            
            # update oligo_conc
            if self.best_bp_distance == 0:
                if self.oligo_conc > target_oligo_conc:
                    if self.oligo_conc/10 <= target_oligo_conc:
                        self.oligo_conc = target_oligo_conc
                    else:
                        self.oligo_conc /= 10
                self.update_sequence(self.sequence)
                self.update_best()
                print self.native, self.bp_distance
                print self.best_native, self.best_bp_distance
            #    self.all_solutions.append([mut_sequence, score])
        
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

