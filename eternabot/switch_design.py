import inv_utils
import eterna_utils
import sys

def check_secstructs(sequence, target, oligo_sequence):
    """
    args:
    sequence is the current sequence of the RNA
    target is dictionary with target secondary structures
    oligo_sequence is the sequence of the oligo to be tested

    return:
    boolean indicating whether the RNA folds to the targetted structure
    with and without the oligo
    """
    secstruct = {}
    secstruct['oligo'] = inv_utils.fold(sequence)
    secstruct['no_oligo'] = inv_utils.fold('&'.join([sequence oligo_sequence]))
    if eterna_utils.bp_distance(secstruct['oligo'], target['oligo']) == 0 and \
        eterna_utils.bp_distance(secstruct['no_oligo'], target['no_oligo']) == 0):
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

    # if secstruct is all unpaired, then keep start sequence    
    dotonly = True
    for ii in range(0,len(secstruct)):
        if(secstruct[ii] != "."):
            dotonly = False
    
    if(dotonly or len(secstruct) < 20): 
         return [start_sequence, 0, 0, {}]  
    
    bases = "GAUC"
    pairs = ["GC", "CG", "AU", "UA"]

    #print "%s %s" % (secstruct, start_sequence)

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
# parse inputs
secstruct = {}
secstruct['oligo'] = sys.argv[1]
secstruct['no_oligo'] = sys.argv[2]
constraints = sys.argv[3]
oligo_sequence = sys.argv[4]
