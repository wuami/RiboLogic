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

# parse inputs
secstruct = {}
secstruct['oligo'] = sys.argv[1]
secstruct['no_oligo'] = sys.argv[2]
constraints = sys.argv[3]
oligo_sequence = sys.argv[4]
