import inv_utils
import eterna_utils
import unittest
import sys

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

def check_secstructs(sequence, target, oligo_sequence):
    """
    args:
    sequence is the current sequence of the RNA
    target is dictionary with target secondary structures and constraint string
    oligo_sequence is the sequence of the oligo to be tested

    return:
    boolean indicating whether the RNA folds to the targetted structure
    with and without the oligo
    """
    l = len(sequence)
    secstruct = {}
    secstruct['oligo'] = inv_utils.cofold('&'.join([sequence, oligo_sequence]))[0][:l]
    secstruct['no_oligo'] = inv_utils.fold(sequence)[0]
    if bp_distance_with_constraint(secstruct['oligo'], target['oligo'][0], target['oligo'][1]) == 0 and \
        bp_distance_with_constraint(secstruct['no_oligo'], target['no_oligo'][0], target['no_oligo'][1]) == 0:
        return True
    return False

def optimize_sequence(secstruct, start_sequence, constraints, scoring_func, score_cutoff):
    """
    args:
    secstruct contains dict where keys are lists of seq numbers and values are
    start_sequence is the current sequence
    constraints contain constraints on sequence positions
    scoring_func is a function for scoring sequences
    score_cutoff is the cutoff score to stop the search
    
    return:
    best sequence, bp distance, final design 
    """

    bases = "GAUC"
    pairs = ["GC", "CG", "AU", "UA"]

    target_pairmap = eterna_utils.get_pairmap_from_secstruct(secstruct)
    n = len(start_sequence)

    # get indices for stacks and loops
    stack_indices = []
    loop_indices = []
    for ii in range(0,len(target_pairmap)):
        if target_pairmap[ii] > ii and constraints[ii] == "N":
            stack_indices.append(ii)
        elif target_pairmap[ii] < 0 and constraints[ii] == "N":
            if(ii > 0 and target_pairmap[ii-1] >= 0):
                loop_indices.append(ii)
            elif(ii < len(target_pairmap)-1 and target_pairmap[ii+1] >= 0):
                loop_indices.append(ii)
    
    # get information about current sequence    
    sequence = start_sequence
    native = inv_utils.fold(sequence)[0]
    bp_distance = eterna_utils.bp_distance(secstruct,native)
    native_pairmap = eterna_utils.get_pairmap_from_secstruct(native)
    design = eterna_utils.get_design_from_sequence(sequence,secstruct)
    design_score = scoring_func(design)
    
    # initialize variables for iteration
    walk_iter = 0
    stale_move = 0

    best_sequence = sequence
    best_bp_distance = bp_distance
    best_native = native
    best_native_pairmap = native_pairmap
    best_design = design
    best_design_score = design_score
    
    # get indices for positions that can change
    index_array = []
    for ii in range(0,n):
        if(constraints[ii] == "N"):
            index_array.append(ii)
    
    # loop as long as bp distance too large or design score too small
    while(bp_distance > 10 or design_score['finalscore'] < score_cutoff) and len(index_array) > 0:
        #print "%d %d %f" %(walk_iter, best_bp_distance, best_design_score['finalscore'])
        #random.shuffle(index_array)
        
        # iterate over nucleotides in sequence
        moved = False
        for rrr in range(0,len(index_array)):
            ii = index_array[rrr]
            
            break_for = False
            if(target_pairmap[ii] != native_pairmap[ii]):
                # nonmatching nucleotides that should be unpaired
                if(target_pairmap[ii] < 0):
                    # try each other possible base in this position
                    for jj in range(0,len(bases)):
                        if sequence[ii] == bases[jj]:
                            continue
                        mut_array = get_sequence_array(sequence)
                        mut_array[ii] = bases[jj]
                        
                        mut_sequence = get_sequence_string(mut_array)
                        mut_native = inv_utils.fold(mut_sequence)[0]
                        mut_bp_distance = eterna_utils.bp_distance(secstruct,mut_native)
                        mut_design = eterna_utils.get_design_from_sequence(mut_sequence,secstruct)
                        mut_score = scoring_func(design)
        
                        # if distance or score is better for mutant, update the current sequence
                        if((mut_bp_distance < bp_distance) or mut_score['finalscore'] > design_score['finalscore']):
                            sequence = mut_sequence
                            native = mut_native
                            bp_distance = mut_bp_distance
                            native_pairmap = eterna_utils.get_pairmap_from_secstruct(native)
                            design = mut_design
                            design_score = mut_score
                        
                            # if distance or score is better for mutant than best, update the best sequence    
                            if((mut_bp_distance < best_bp_distance or mut_bp_distance < 10 or len(secstruct) < 50) and mut_score['finalscore'] > best_design_score['finalscore']):
                                best_sequence = mut_sequence
                                best_bp_distance = mut_bp_distance
                                best_native = mut_native
                                best_native_pairmap = eterna_utils.get_pairmap_from_secstruct(mut_native)
                                best_design = mut_design
                                best_design_score = mut_score
                            
                            break_for = True
                            break
                # nonmatching nucleotides that should be paired
                else:
                    # try each other possible pair in this position
                    current_pair = sequence[ii] + sequence[target_pairmap[ii]]
                    for jj in range(0,len(pairs)):
                        if current_pair == pairs[jj]:
                            continue
                        mut_array = get_sequence_array(sequence)
                        mut_array[ii] = pairs[jj][0]
                        mut_array[target_pairmap[ii]] = pairs[jj][1]
                        
                        mut_sequence = get_sequence_string(mut_array)
                        mut_native = inv_utils.fold(mut_sequence)[0]
                        mut_bp_distance = eterna_utils.bp_distance(secstruct,mut_native)
                        mut_design = eterna_utils.get_design_from_sequence(mut_sequence,secstruct)
                        mut_score = scoring_func(design)
                        
                        # if distance or score is better for mutant, update the current sequence
                        if((mut_bp_distance < bp_distance) or mut_score['finalscore'] > design_score['finalscore']):
                            sequence = mut_sequence
                            bp_distance = mut_bp_distance
                            native = mut_native
                            native_pairmap = eterna_utils.get_pairmap_from_secstruct(native)
                            design = mut_design
                            design_score = mut_score
                            
                            # if distance or score is better for mutant than best, update best
                            if((mut_bp_distance < best_bp_distance or mut_bp_distance < 10 or len(secstruct) < 50) and mut_score['finalscore'] > best_design_score['finalscore']):
                                best_sequence = mut_sequence
                                best_bp_distance = mut_bp_distance
                                best_native = mut_native
                                best_native_pairmap = eterna_utils.get_pairmap_from_secstruct(mut_native)
                                best_design = mut_design
                                best_design_score = mut_score
                            
                            break_for = True
                            break
            # if the sequence has updated, break loop through nucleotides   
            if(break_for):
                moved = True
                break
        
        # if sequence hasn't been updated, randomize sequence
        if(moved == False):

            stale_move += 1
            
            if(stale_move > 10):
                return [best_sequence, best_bp_distance, best_design_score['finalscore'], best_design_score]
            
            rand = random.random()
        
            mut_array = get_sequence_array(sequence)    
            if(rand < 0.5 and len(stack_indices) > 0):              
                rindex = stack_indices[int(random.random() * len(stack_indices)) % len(stack_indices)]
                if(target_pairmap[rindex] < 0):
                    print "Something is wrong"
                    sys.exit(0)
                
                if(rand < 0.33):
                    temp = mut_array[rindex]
                    mut_array[rindex] = mut_array[target_pairmap[rindex]]
                    mut_array[target_pairmap[rindex]] = temp
                else:
                    pair = get_random_pair()
                    mut_array[rindex] = pair[0]
                    mut_array[target_pairmap[rindex]] = pair[1]
            else:
                if (len(loop_indices) > 0):
                    rindex = loop_indices[int(random.random() * len(loop_indices)) % len(loop_indices)]
                    mut_array[rindex] = get_random_base()
            
            sequence = get_sequence_string(mut_array)
            native = inv_utils.fold(sequence)[0]
            bp_distance = eterna_utils.bp_distance(secstruct,native)
            design = eterna_utils.get_design_from_sequence(sequence,secstruct)
            design_score = scoring_func(design)
        else:
            stale_move = 0
            
        moved == False
        walk_iter += 1
    
        # if it has been too many iterations, finish    
        if(walk_iter > 400):
            return [best_sequence, best_bp_distance, best_design_score['finalscore'], best_design_score]
    
    return [best_sequence, best_bp_distance, best_design_score['finalscore'], best_design_score]                

class test_functions(unittest.TestCase):

    def setUp(self):
        with open('switch_input.txt', 'r') as f:
            self.secstruct = {}
            self.secstruct['oligo'] = f.readline().split()
            self.secstruct['no_oligo'] = f.readline().split()
            self.constraints = f.readline().strip()
            self.oligo_sequence = f.readline().strip()

    def test_check_secstructs(self):
        sequence = "ACAAGCUUUUUGCUCGUCUUAUACAUGGGUAAAAAAAAAACAUGAGGAUCACCCAUGUAAAAAAAAAAAAAAAAAAA"
        self.assertTrue(check_secstructs(sequence, self.secstruct, self.oligo_sequence))

if __name__ == "__main__":
    unittest.main()
