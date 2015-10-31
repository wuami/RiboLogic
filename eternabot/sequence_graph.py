import ensemble_design, design_utils, eterna_utils
import random
import networkx as nx
import matplotlib.pyplot as plt
import sys
import itertools

base_coloring = {'A':'y', 'U':'b', 'G':'r', 'C':'g', '&':'w'}

class SequenceGraph(object):

    def __init__(self, inputs, targets, seq_locks, sequence, add_rcs, draw=False, autocomplement=True):
        """
        create graph structure to represent sequence
        given inputs, target secondary structures, sequence constraints
        """

        # set class variables
        self.inputs = inputs
        self.targets = targets
        self.seq_locks = self._get_full_seqlocks(seq_locks)
        self.design_seq = sequence
        self.sequence = self._get_full_sequence(sequence)
        self.n = len(sequence)
        self.N = len(self.sequence)
        self.autocomplement = autocomplement

        # create dependency graph
        self.index_array = self.get_unconstrained_indices()
        self.dep_graph = self.get_dependency_graph()
        #print self.dep_graph.nodes(data=True)

        # update constrained parts of sequence
        self.add_rcs = add_rcs
        self.reset_sequence(self.design_seq)

        # draw if option specified
        self.draw = draw
        if self.draw:
            self.printi = 0
            # draw dependency graph, for debugging
            nx.draw_circular(self.dep_graph, with_labels=True, node_color=[base_coloring[base] for base in self.sequence], node_size=60, font_size=5)
            plt.savefig("dependency_graph%s.png" % self.printi, dpi=300)
            self.printi += 1
        
        # set mutation function
        if add_rcs:
            self.mutate_func = self.mutate_or_shift
        else:
            self.mutate_func = self.mutate_sequence

    def reset_sequence(self, sequence):
        """
        set sequence back to inital
        """
        self.sequence = self._get_full_sequence(sequence)
        seq_array = ensemble_design.get_sequence_array(self.sequence)
        if self.autocomplement:
            for i in range(self.N):
                if self.seq_locks[i] == "x":
                    self.update_neighbors(i, seq_array, [])
        self.sequence = ensemble_design.get_sequence_string(seq_array)
        if self.add_rcs:
            self.init_oligo_rcs()

    def init_oligo_rcs(self):
        """
        set reverse complements to inputs
        """
        self.oligo_rc = []
        self.oligo_pos = []
        self.oligo_len = []
        self.oligo_len_sum = 0
        inputs = self.inputs.values()
        oligo1_start = [x + self.seq_offset for x in [0, len(inputs[0])]]
        self.set_oligo_rc(inputs[0], oligo1_start)
        oligo2_start = [x + self.seq_offset for x in [self.n-len(inputs[1]), self.n]]
        self.set_oligo_rc(inputs[1], oligo2_start)
        self.pmut = 0.2
        self.pexpand = 0.8
        self.pshift = 0.5

    def set_oligo_rc(self, seq, range):
        """
        set reverse complement parameters for one input
        """
        #lo, hi = self.n-49, self.n-27
        rc = design_utils.rc(seq)
        self.sequence = self.sequence[0:range[0]] + rc + self.sequence[range[1]:]
        self.oligo_rc.append(rc)
        self.oligo_pos.append(range)
        self.oligo_len.append([0,range[1]-range[0]-1])
        self.oligo_len_sum += len(seq)

    def _get_full_secstruct(self, secstruct, inputs):
        """
        get full secondary structure string with inputs
        """
        substrings = secstruct.split('&')
        secstruct = ""
        i = 0
        for input in sorted(self.inputs):
            if input in inputs:
                secstruct += substrings[i] + '&'
                i += 1
            else:
                secstruct += '.'*len(self.inputs[input]) + '&'
        return secstruct + substrings[-1]
        
    def _get_full_seqlocks(self, seq_locks):
        """
        get full sequence constraint string with inputs
        """
        constraint = ""
        for input in sorted(self.inputs):
            constraint += "x"*len(self.inputs[input]) + "o"
        return constraint + seq_locks
    
    def _get_full_sequence(self, sequence):
        """
        get full sequence string with inputs
        """
        seq = ""
        for input in sorted(self.inputs):
            seq += self.inputs[input] + "&"
        return seq + sequence

    def get_dependency_graph(self):
        """ 
        get dependency graph based on target secondary structures
        """
        # position of start of actual sequence in graph
        self.seq_offset = sum([len(x) for x in self.inputs.values()]) + len(self.inputs)

        # create graph
        graph = nx.Graph()
        graph.add_nodes_from(range(self.N), bases=['A','U','G','C'])
        for target in self.targets:
            if "inputs" in target:
                target['full_secstruct'] = self._get_full_secstruct(target['secstruct'], target['inputs'])
            else:
                target['full_secstruct'] = target['secstruct']
            pairmap = eterna_utils.get_pairmap_from_secstruct(target['full_secstruct'])
            for i, j in enumerate(pairmap):
                if j > i:
                    graph.add_edge(i,j)
        if not nx.is_bipartite(graph):
            raise ValueError("dependency graph is not bipartite")
        self.dep_graph = graph
        if self.autocomplement:
            self.set_sequence_constraints()
        return graph

    def set_sequence_constraints(self):
        """
        set sequence constraints for all positions based on graph structure
        """
        for i in range(len(self.seq_locks)):
            if self.seq_locks[i] == "x":
                self.set_neighbor_sequence_constraints(i, [self.sequence[i]], [])
    
    def set_neighbor_sequence_constraints(self, i, bases, updated):
        """
        recursively set sequence constraints for position i
        """
        self.dep_graph.node[i]['bases'] = list(set(self.dep_graph.node[i]['bases']) & set(bases))
        updated.append(i)
        for pos in nx.all_neighbors(self.dep_graph, i):
            if pos not in updated:
                possible_bases = []
                for base in bases:
                    possible_bases.extend(design_utils.get_rcs(base))
                self.set_neighbor_sequence_constraints(pos, possible_bases, updated)
            
    def get_unconstrained_indices(self):
        """
        get indices for positions that can change
        """
        open = []
        restricted = []
        for ii in range(0,self.N):
            if(self.seq_locks[ii] == "o"):
                open.append(ii)
            elif self.seq_locks[ii] == "r":
                restricted.append(ii)
        return [open, restricted]

    def mutate(self):
        """
        general mutate function
        """
        self.sequence = self.mutate_func()
        if self.draw:
            if self.printi <= 20:
                ## draw dependency graph, for debugging
                nx.draw_circular(self.dep_graph, with_labels=True, node_color=[base_coloring[base] for base in self.sequence], node_size=60, font_size=5)
                plt.savefig("dependency_graph%s.png" % self.printi, dpi=300)
                self.printi += 1
        return self.sequence[-self.n:]

    def get_rand_shift(self):
        # randomly choose an oligo and whether to shift
        roligo = random.getrandbits(1)
        shift = random.random() < self.pshift
        min_index = min(list(itertools.chain(*self.index_array)))
        max_index = max(list(itertools.chain(*self.index_array)))

        # choose to expand or shrink rc sequence
        maxed_out = False
        if shift:
            expand = 0
        else:
            if self.oligo_len[roligo][1]-self.oligo_len[roligo][0] <= 0:
                expand = 1
            if self.oligo_len[roligo][1]-self.oligo_len[roligo][0] == len(self.oligo_rc)-1:
                expand = 0
            if 'expand' not in locals():
                expand = random.random() < self.pexpand


        # choose left or right
        if expand or shift:
            if self.oligo_pos[roligo][1] >= self.N or self.oligo_pos[roligo][1] >= max_index or \
               (self.oligo_len[roligo][1] >= len(self.oligo_rc[roligo])-1 and expand) or \
               self.seq_locks[self.oligo_pos[roligo][1]-1] == 'x':
                right = 0
                maxed_out = True
            if self.oligo_pos[roligo][0] <= 0 or self.oligo_pos[roligo][0] <= min_index or \
               (self.oligo_len[roligo][0] <= len(self.oligo_rc[roligo])/2 and expand) or \
               self.seq_locks[self.oligo_pos[roligo][0]-1] == 'x':
                right = 1
                if maxed_out:
                    expand = 0
        if 'right' not in locals():
            right = random.getrandbits(1)
        
        return roligo, shift, expand, right

    def mutate_or_shift(self):
        """
        either mutate sequence or shift the complement to the oligo in the design sequence
        """
        # mutate randomly wp 0.5, otherwise mutate oligo rc
        if (random.random() < self.pmut): #float(self.oligo_len_sum)/sum([len(x) for x in self.index_array])):
            return self.mutate_sequence()

        mut_array = ensemble_design.get_sequence_array(self.sequence)
        roligo, shift, expand, right = self.get_rand_shift()
        # shift complement sequence
        if shift:
            if right:
                mut_array[self.oligo_pos[roligo][0]] = ensemble_design.get_random_base(self.dep_graph.node[self.oligo_pos[roligo][0]]['bases'])
                self.oligo_pos[roligo] = [x+1 for x in self.oligo_pos[roligo]] 
            else:
                self.oligo_pos[roligo] = [x-1 for x in self.oligo_pos[roligo]] 
                mut_array[self.oligo_pos[roligo][1]] = ensemble_design.get_random_base(self.dep_graph.node[self.oligo_pos[roligo][1]]['bases'])
            seq = ensemble_design.get_sequence_string(mut_array)
            start, end = self.oligo_pos[roligo]
            lo, hi = self.oligo_len[roligo]
            new_seq = seq[:start] + self.oligo_rc[roligo][lo:hi+1] + seq[end:]
            assert len(seq) == len(new_seq), "shifted sequence changed lengths"
            return new_seq
        # expand complement sequence
        if expand:
            if right:
                mut_array[self.oligo_pos[roligo][right]] = self.oligo_rc[roligo][self.oligo_len[roligo][right]]
                self.oligo_pos[roligo][1] += 1
                self.oligo_len[roligo][1] += 1
            else:
                self.oligo_pos[roligo][0] -= 1
                self.oligo_len[roligo][0] -= 1
                mut_array[self.oligo_pos[roligo][right]] = self.oligo_rc[roligo][self.oligo_len[roligo][right]]
        # shrink complement sequence
        else:
            if right:
                self.oligo_pos[roligo][1] -= 1
                self.oligo_len[roligo][1] -= 1
                mutate_pos = self.oligo_pos[roligo][right]
                mut_array[mutate_pos] = ensemble_design.get_random_base(self.dep_graph.node[mutate_pos]['bases'])
            else:
                mutate_pos = self.oligo_pos[roligo][right]
                mut_array[mutate_pos] = ensemble_design.get_random_base(self.dep_graph.node[mutate_pos]['bases'])
                self.oligo_pos[roligo][0] += 1
                self.oligo_len[roligo][0] += 1 
        self.oligo_len_sum = sum([x[1]-x[0] for x in self.oligo_len])
        print self.oligo_pos, self.oligo_len
        return ensemble_design.get_sequence_string(mut_array)

    def mutate_sequence(self):
        """ mutate one random position """
        while True:
            mut_array = ensemble_design.get_sequence_array(self.sequence)
            if random.random() < 0.9 or len(self.index_array[1]) == 0 and len(self.index_array[0]) != 0:
                rindex = self.index_array[0][int(random.random() * len(self.index_array[0]))]
            else:
                rindex = self.index_array[1][int(random.random() * len(self.index_array[1]))]
            rbase = ensemble_design.get_random_base(self.dep_graph.node[rindex]['bases'])
            mut_array[rindex] = rbase
            if self.autocomplement:
                self.update_neighbors(rindex, mut_array, [])
            if design_utils.satisfies_constraints(mut_array, self.sequence, self.seq_locks):
                break
        return ensemble_design.get_sequence_string(mut_array)

    def update_neighbors(self, node, mut_array, updated):
        """
        recursively update neighbors to specified node
        """
        for pos in nx.all_neighbors(self.dep_graph, node):
            if pos not in updated:
                complement = design_utils.rc(mut_array[node], possible_bases=self.dep_graph.node[pos]['bases'])
                mut_array[pos] = complement
                updated.append(pos)
                self.update_neighbors(pos, mut_array, updated)
