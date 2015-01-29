import ensemble_design
import ensemble_utils
import eterna_utils
import switch_designer, gate_designer
import sys
import os
import json
import argparse
import requests
import settings

def read_puzzle_json(text):
    """
    read in puzzle as a json file
    """
    data = json.loads(text)['data']
    p = data['puzzle']
    id = data['nid']

    #if p['rna_type'] == 'single':
    #    return SequenceDesigner(id, p['secstruct'], p['locks'])

    # get basic parameters
    beginseq = p['beginseq']
    constraints = p['locks']

    # load in objective secondary structures
    objective = json.loads(p['objective'])
    secstruct = [] 
    for o in objective:
        n = len(o['secstruct'])
        # if no constrained bases, all are unconstrained
        if 'structure_constrained_bases' not in o.keys() and 'anti_structure_constrained_bases' not in o.keys():
            constrained = ensemble_design.get_sequence_array('x'*n)
        # otherwise, combine structure and antistructure constraints
        else:
            constrained = ensemble_design.get_sequence_array('o'*n)
            struct = ensemble_design.get_sequence_array(o['secstruct'])
            if 'structure_constrained_bases' in o.keys() and len(o['structure_constrained_bases']) > 0:
                [lo, hi] = o['structure_constrained_bases']
                for i in range(lo, hi+1):
                    constrained[i] = 'x'
                del o['structure_constrained_bases']
            if 'anti_structure_constrained_bases' in o.keys() and len(o['anti_structure_constrained_bases']) > 0:
                [lo, hi] = o['anti_structure_constrained_bases']
                for i in range(lo, hi+1):
                    constrained[i] = 'x'
                    struct[i] = '.'
                del o['anti_secstruct'], o['anti_structure_constrained_bases']
            o['secstruct'] = ensemble_design.get_sequence_string(struct)
        o['constrained'] = ensemble_design.get_sequence_string(constrained)
        secstruct.append(o)

    # create scoring function
    strategy_names = ['merryskies_only_as_in_the_loops', 'aldo_repetition', 'dejerpha_basic_test', 'eli_blue_line', 'clollin_gs_in_place', 'quasispecies_test_by_region_boundaries', 'eli_gc_pairs_in_junction', 'eli_no_blue_nucleotides_in_hook', 'mat747_31_loops', 'merryskies_1_1_loop', 'xmbrst_clear_plot_stack_caps_and_safe_gc', 'jerryp70_jp_stratmark', 'eli_energy_limit_in_tetraloops', 'eli_double_AUPair_strategy', 'eli_green_blue_strong_middle_half', 'eli_loop_pattern_for_small_multiloops', 'eli_tetraloop_similarity', 'example_gc60', 'penguian_clean_dotplot', 'eli_twisted_basepairs', 'aldo_loops_and_stacks', 'eli_direction_of_gc_pairs_in_multiloops_neckarea', 'eli_multiloop_similarity', 'eli_green_line', 'ding_quad_energy', 'quasispecies_test_by_region_loops', 'berex_berex_loop_basic', 'eli_legal_placement_of_GUpairs', 'merryskies_1_1_loop_energy', 'ding_tetraloop_pattern', 'aldo_mismatch', 'eli_tetraloop_blues', 'eli_red_line', 'eli_wrong_direction_of_gc_pairs_in_multiloops', 'deivad_deivad_strategy', 'eli_direction_of_gc_pairs_in_multiloops', 'eli_no_blue_nucleotides_strategy', 'berex_basic_test', 'eli_numbers_of_yellow_nucleotides_pr_length_of_string', 'kkohli_test_by_kkohli']
    weights_file_name = "no_validation_training/weights_sparse_5.overall.txt"
    scores_file_name = "no_validation_training/predicted_score_sparse_5.overall.unnormalized.txt"
    weights_f = open(os.path.join(settings.RESOURCE_DIR, weights_file_name),"r")
    weights = []
    for line in weights_f:
        weights.append(float(line))
    ensemble = ensemble_utils.Ensemble("sparse", strategy_names, weights)

    if p['rna_type'] == "multi_input":
        return switch_designer.SwitchDesigner(id, p['rna_type'], beginseq, constraints, secstruct, ensemble.score, p['inputs'])
    elif p['rna_type'] == "multi_input_oligo":
        return gate_designer.GateDesigner(id, p['rna_type'], beginseq, constraints, secstruct, ensemble.score, p['inputs'])
    return switch_designer.SwitchDesigner(id, p['rna_type'], beginseq, constraints, secstruct, ensemble.score)
    
def optimize_n(puzzle, niter, ncool, n, submit, fout):
    if fout:
        with open(fout, 'a') as f:
	        f.write("# %s iterations, %s coolings\n" % (niter, ncool))

    # run puzzle n times
    solutions = []
    scores = []
    i = 0 
    attempts = 0
    while i < n:
        puzzle.reset_sequence()
        puzzle.optimize_sequence(niter, ncool)
        if puzzle.check_current_secstructs():
            sol = puzzle.get_solution()
            if sol[0] not in solutions:
                solutions.append(sol[0])
                scores.append(sol[2])
                print sol
                if submit:
                    post_solution(puzzle, 'solution %s' % i)
                if fout:
                    with open(fout, 'a') as f:
                        f.write("%s\t%1.6f\n" % (sol[0], sol[2]))
                i += 1
                attempts = 0
        else:
            #niter += 500
            attempts += 1
            if attempts == 10:
                break
        print "%s sequence(s) calculated" % i
    return [solutions, scores]

def get_puzzle(id):
    """
    get puzzle with id number id from eterna server
    """
    r = requests.get('http://nando.eternadev.org/get/?type=puzzle&nid=%s' % id)
    return read_puzzle_json(r.text)

def post_solution(puzzle, title):
    sequence = puzzle.best_sequence
    fold = inv_utils.fold(sequence)
    design = eterna_utils.get_design_from_sequence(sequence, fold[0])
    header = {'Content-Type': 'application/x-www-form-urlencoded'}
    login = {'type': 'login',
             'name': 'theeternabot',
             'pass': 'iamarobot',
             'workbranch': 'main'}
    solution = {'type': 'post_solution',
                'puznid': puzzle.id,
                'title': title,
                'body': 'eternabot switch v1, score %s' % puzzle.best_design_score,
                'sequence': sequence,
                'energy': fold[1],
                'gc': design['gc'],
                'gu': design['gu'],
                'ua': design['ua'],
                'melt': design['meltpoint'],
                'pointsrank': 'false',
                'recommend-puzzle': 'true'}

    url = "http://jnicol.eternadev.org"
    #url = 'http://eterna.cmu.edu'
    loginurl = "%s/login/" % url
    posturl = "%s/post/" % url
    with requests.Session() as s:
        r = s.post(loginurl, data=login, headers=header)
        r = s.post(posturl, data=solution, headers=header)
    return

def main():
    p = argparse.ArgumentParser()
    p.add_argument('puzzleid', help="name of puzzle filename or eterna id number", type=str)
    p.add_argument('-s', '--nsol', help="number of solutions", type=int, default=1)
    p.add_argument('-i', '--niter', help="number of iterations", type=int, default=1000)
    p.add_argument('-c', '--ncool', help="number of times to cool", type=int, default=20)
    p.add_argument('--submit', help="submit the solution(s)", default=False, action='store_true')
    p.add_argument('--nowrite', help="suppress write to file", default=False, action='store_true')
    args = p.parse_args()

    puzzlefile = os.path.join(settings.PUZZLE_DIR, "%s.json" % args.puzzleid)
    if os.path.isfile(puzzlefile): 
        with open(puzzlefile, 'r') as f:
            puzzle = read_puzzle_json(f.read())
    else:
        puzzle = get_puzzle(args.puzzleid)
    if not args.nowrite:
        fout = os.path.join(settings.PUZZLE_DIR, args.puzzleid + ".out")
    else:
        fout = False
    [solutions, scores] = optimize_n(puzzle, args.niter, args.ncool, args.nsol, args.submit, fout)

if __name__ == "__main__":
    #unittest.main()
    main()


