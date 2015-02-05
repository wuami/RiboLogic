import eterna_utils
import inv_utils
import random
import math
import ensemble_design
import switch_designer

class GateDesigner(switch_designer.SwitchDesigner):

    def __init__(self, id, type, beginseq, constraints, targets, scoring_func, inputs = None):
        switch_designer.SwitchDesigner.__init__(self, id, type, beginseq, constraints, targets, scoring_func, inputs)
        self.set_oligo_rc()

    def set_oligo_rc(self):
        self.oligo_rc = self.beginseq[27:49]
        self.oligo_pos = [27,48]
        #self.oligo_rc = self.beginseq[55:77]
        #self.oligo_pos = [55,76]
        self.oligo_len = [0,21]

    def reset_sequence(self):
        super(GateDesigner, self).reset_sequence()
        self.set_oligo_rc()

    def score_secstructs(self, secstruct):
        """
        calculates sum of bp distances for with and without oligo

        returns:
        sum of bp distances with and without the oligo 
        """
        distance = 0
        for i in range(self.n_targets):
            if "threshold" in self.targets[i]:
                dist = self.bp_distance_func(secstruct[i], self.targets[i]['secstruct'], self.targets[i]['constrained'], self.targets[i]['threshold'])
            else:
                dist = self.bp_distance_func(secstruct[i], self.targets[i]['secstruct'], self.targets[i]['constrained'])
            distance += dist
        return distance

    def mutate_sequence(self, sequence):
        """
        either mutate sequence or shift the complement to the oligo in the design sequence
        """
        mut_array = ensemble_design.get_sequence_array(self.sequence)
        # mutate randomly wp 0.5, otherwise mutate oligo rc
        if (random.random() > float(self.oligo_len[1]-self.oligo_len[0]+1)/len(self.index_array)):
            rindex = self.index_array[int(random.random() * len(self.index_array))]
            #while rindex in range(self.oligo_pos[0], self.oligo_pos[1]+1):
            #    rindex = self.index_array[int(random.random() * len(self.index_array))]
            mut_array[rindex] = ensemble_design.get_random_base()
        else:
            ## wp 0.5 change length
            #if random.random() < 0.5:
            rindex = random.getrandbits(1) # pick left or right
            # wp 0.5 expand
            if (random.random() < 0.5 or self.oligo_len[1]-self.oligo_len[0] <= 0) and \
                self.oligo_len[1]-self.oligo_len[0] != len(self.oligo_rc)-1:
                if (rindex or self.oligo_pos[0] == self.index_array[0] or self.oligo_len[0] == 0) and \
                    self.oligo_pos[1] != self.index_array[-1] and self.oligo_len[1] != len(self.oligo_rc)-1 and\
                    self.constraints[self.oligo_pos[rindex]+1] != 'x':
                    self.oligo_pos[rindex] += 1
                    self.oligo_len[rindex] += 1
                elif self.constraints[self.oligo_pos[rindex]-1] != 'x':
                    self.oligo_pos[rindex] -= 1
                    self.oligo_len[rindex] -= 1
                mut_array[self.oligo_pos[rindex]] = self.oligo_rc[self.oligo_len[rindex]]
            # otherwise shrink
            else:
                mut_array[self.oligo_pos[rindex]] = ensemble_design.get_random_base()
                if rindex:
                    self.oligo_pos[rindex] -= 1
                    self.oligo_len[rindex] -= 1
                else:
                    self.oligo_pos[rindex] += 1
                    self.oligo_len[rindex] += 1 
            ## otherwise shift
            #else:
            #    # wp 0.5 shift left
            #    if (random.random() < 0.5 or self.oligo_pos[1] == self.n-1) and self.oligo_pos[0] != 0:
            #        print "shift left"
            #        self.oligo_pos = [i-1 for i in self.oligo_pos]
            #        for i in range(self.oligo_len[0], self.oligo_len[1]+1):
            #            mut_array[self.oligo_pos[0]+i] = self.oligo_rc[i]
            #        mut_array[self.oligo_pos[1]+1] = ensemble_design.get_random_base()
            #    # otherwise shift right
            #    else:
            #        print "shift right"
            #        self.oligo_pos = [i+1 for i in self.oligo_pos]
            #        for i in range(self.oligo_len[0], self.oligo_len[1]+1):
            #            mut_array[self.oligo_pos[0]+i] = self.oligo_rc[i]
            #        mut_array[self.oligo_pos[0]-1] = ensemble_design.get_random_base()
        return ensemble_design.get_sequence_string(mut_array)


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

        #for target in self.targets:
        #    print target

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
            mut_sequence = self.mutate_sequence(self.sequence)
            [mut_sequence, native, bp_distance, score] = self.get_sequence_info(mut_sequence)

            # if current sequence is a solution, save to list
            if bp_distance == 0:
                self.all_solutions.append([mut_sequence, score])
            
            # if distance or score is better for mutant, update the current sequence
            if(random.random() < p_dist(self.bp_distance, bp_distance) or
               (bp_distance == self.bp_distance and random.random() < p_score(self.design_score, score))):
                self.update_sequence(mut_sequence, native, bp_distance, score)
                #print self.sequence, self.bp_distance
                #for struct in native:
                #    print struct[0:22], struct[78:97]
            
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

