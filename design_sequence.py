import switch_designer
import sys, os
import json, re
import argparse
import requests
import settings
import varna
import copy
import signal

def get_objective_dict(o):
    n = len(o['secstruct'])
    constrained = ['o']*n
    struct = list(o['secstruct'])
    if 'anti_structure_constrained_bases' in o.keys() and len(o['anti_structure_constrained_bases']) > 0:
        for i in range(0, len(o['anti_structure_constrained_bases']), 2):
            [lo, hi] = o['anti_structure_constrained_bases'][i:i+2]
            for j in range(lo, hi+1):
                constrained[j] = 'x'
                struct[j] = '.'
        del o['anti_secstruct'], o['anti_structure_constrained_bases']
    if 'structure_constrained_bases' in o.keys() and len(o['structure_constrained_bases']) > 0:
        for i in range(0, len(o['structure_constrained_bases']), 2):
            [lo, hi] = o['structure_constrained_bases'][i:i+2]
            for j in range(lo, hi+1):
                constrained[j] = 'x'
                struct[j] = o['secstruct'][j]
        del o['structure_constrained_bases']
    if 'structure_unpaired_constrained_bases' in o.keys() and len(o['structure_unpaired_constrained_bases']) > 0:
        for i in range(0, len(o['structure_unpaired_constrained_bases']), 2):
            [lo, hi] = o['structure_unpaired_constrained_bases'][i:i+2]
            for j in range(lo, hi+1):
                constrained[j] = 'u'
        del o['structure_unpaired_constrained_bases']
    if 'structure_paired_constrained_bases' in o.keys() and len(o['structure_paired_constrained_bases']) > 0:
        for i in range(0, len(o['structure_paired_constrained_bases']), 2):
            [lo, hi] = o['structure_paired_constrained_bases'][i:i+2]
            for j in range(lo, hi+1):
                constrained[j] = 'p'
        del o['structure_paired_constrained_bases']
    o['secstruct'] = ''.join(struct)
    o['constrained'] = ''.join(constrained)
    return o

def read_constraints_from_file(filename, **kwargs):
    """
    reads design information from text file
    """
    inputs = {}
    targets = []
    beginseq = ''
    constraints = ''
    
    with open(filename) as f:
        line = f.readline()
        while line:
            # read inputs
            if line.startswith('<'):
                input = f.readline().strip()
                try:
                    inputs[line.strip('<\n')] = {'type':'ligand', 'kD': float(input), 'fold_constraint': f.readline().strip()}
                except:
                    inputs[line.strip('<\n')] = {'type':'RNA', 'sequence':f.readline().strip()}
            # read sequence constraints
            elif line.startswith('-'):
                seq = f.readline().strip()
                beginseq = seq.replace('N', 'A')
                constraints = ''.join(['o' if c == 'N' else 'x' for c in seq])
            # read objectives
            elif line.startswith('>'):
                target = {}
                target['type'] = line.strip('>\n')
                if target['type'] != 'single':
                    target['inputs'] = {}
                    line = f.readline()
                    if not line.strip() == "":
                        for input in line.split(';'):
                            spl = input.split()
                            if len(spl) > 1:
                                target['inputs'][spl[0]] = float(spl[1])
                            else:
                                target['inputs'][spl[0]] = 1
                target['secstruct'] = f.readline().strip()
                target['constrained'] = f.readline().strip()
                line = f.readline()
                if not line.strip() == "":
                    thresholds = [int(x) for x in line.split()]
                    target['threshold'] = []
                    r = re.compile("[up]+[ox]")
                    for i,substruct in enumerate(r.finditer(target['constrained'])):
                        target['threshold'].append([substruct.start(), 
                                                    substruct.start() + len(substruct.group()) - 1,
                                                    int(thresholds[i])])
                targets.append(target)
            line = f.readline()

    return switch_designer.SwitchDesigner(os.path.basename(filename).split('.')[0], beginseq, constraints, targets, inputs=inputs, **kwargs)
    
def read_puzzle_json(text, **kwargs):
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

    outputs = False
    if 'outputs' in p:
        outputs = {}
        for key in p['outputs']:
            outputs[key] = get_objective_dict(p['outputs'][key])

    # load in objective secondary structures
    objective = json.loads(p['objective'])
    secstruct = [] 
    for o in objective:
        if outputs:
            objective = copy.deepcopy(outputs[o['output']])
            objective['type'] = o['type']
            objective['inputs'] = o['inputs']
            secstruct.append(objective)
        else:
            secstruct.append(get_objective_dict(o))

    if 'linker' not in p:
        p['linker'] = 'AACAA'

    if p['rna_type'] == 'multi_input' or p['rna_type'] == 'multi_input_oligo':
        kwargs['inputs'] = p['inputs']
    return switch_designer.SwitchDesigner(id, beginseq, constraints, secstruct, **kwargs)

def optimize_n(puzzle, niter, ncool, n, **kwargs):
    # run puzzle n times
    solutions = []
    scores = []
    i = 0 
    attempts = 0
    while i < n:
        puzzle.reset_sequence()
        passkwargs = {key:kwargs[key] for key in ['greedy', 'start_oligo_conc']}
        nfin = puzzle.optimize_sequence(niter, ncool, **passkwargs)
        if puzzle.check_current_secstructs():
            sol = puzzle.get_solution()
            if sol[0] not in solutions:
                solutions.append(sol[0])
                scores.append(sol[2])
                if 'draw' in kwargs and kwargs['draw']:
                    puzzle.draw_solution(i)
                if 'fout' in kwargs and kwargs['fout']:
                    params = ''
                    if 'greedy' in kwargs and kwargs['greedy']:
                        params += 'greedy '
                    with open(kwargs['fout'], 'a') as f:
                        f.write('# %s out of %s iterations, %s coolings, %s\n' % (nfin, niter, ncool, params))
                        f.write('%s\t%1.6f\n' % (sol[0], sol[2]))
                i += 1
                attempts = 0
        else:
            #niter += 500
            print 'best distance: %s' % puzzle.best_bp_distance
            print 'final conc: %s' % puzzle.oligo_conc
            attempts += 1
            if attempts == 10:
                break
        print '%s sequence(s) calculated' % i
    return [solutions, scores]

def optimize_timed(puzzle, niter, ncool, time, **kwargs):
    def handler(signum, frame):
        raise Exception('%d elapsed' % time)

    # run puzzle n times
    solutions = []
    scores = []
    niters = []
    i = 0 
    signal.signal(signal.SIGALRM, handler)
    signal.alarm(time)
    try:
        while True:
            puzzle.reset_sequence()
            passkwargs = {key:kwargs[key] for key in ['greedy', 'start_oligo_conc']}
            n = puzzle.optimize_sequence(niter, ncool, **passkwargs)
            if puzzle.check_current_secstructs():
                sol = puzzle.get_solution()
                solutions.append(sol[0])
                scores.append(sol[2])
                niters.append(n)
                i += 1
            else:
                print 'best distance: %s' % puzzle.best_bp_distance
                print 'final conc: %s' % puzzle.oligo_conc
            print '%s sequence(s) calculated' % i
    except Exception, exc:
        print exc
        print '%s sequence(s) calculated in %d seconds' % (i, time)
        print 'average iterations: %d' % (float(sum(niters))/len(niters))
        for i in range(len(solutions)):
            print '\t%s %d %d' % (solutions[i], scores[i], niters[i])
    return [solutions, scores]

def get_puzzle(id, **kwargs):
    puzzlefile = os.path.join(settings.PUZZLE_DIR, '%s.json' % id)
    with open(puzzlefile, 'r') as f:
        puzzle = read_puzzle_json(f.read(), **kwargs)
    return puzzle

def view_sequence(puzzle, seq):
    puzzle.update_sequence(seq)
    puzzle.update_best()
    print puzzle.targets
    print puzzle.get_solution()

def main():
    # parse arguments
    p = argparse.ArgumentParser()
    p.add_argument('puzzleid', help='name of puzzle filename or eterna id number', type=str)
    p.add_argument('-n', '--nsol', help='number of solutions', type=int, default=1)
    p.add_argument('-t', '--time', help='maximum time allowed', type=int)
    p.add_argument('-i', '--niter', help='number of iterations', type=int, default=2000)
    p.add_argument('-o', '--ncool', help='number of times to cool', type=int, default=50)
    p.add_argument('-m', '--mode', help='mode for multi inputs', type=str, default='nupack')
    p.add_argument('-c', '--conc', help='starting oligo concentration', type=float, default=1)
    p.add_argument('--draw', help='draw the solution(s)', default=False, action='store_true')
    p.add_argument('--nowrite', help='suppress write to file', default=False, action='store_true')
    p.add_argument('--print_', help='print sequences throughout optimization', default=False, action='store_true')
    p.add_argument('--greedy', help='greedy search', default=False, action='store_true')
    p.add_argument('--add_rcs', help='introduce reverse complement of input oligos', default=False, action='store_true')
    args = p.parse_args()

    print args.puzzleid

    # read puzzle
    #puzzle = get_puzzle(args.puzzleid, mode=args.mode, add_rcs=args.add_rcs, print_=args.print_)
    puzzle = read_constraints_from_file(args.puzzleid, mode=args.mode, add_rcs=args.add_rcs, print_=args.print_)
    if not args.nowrite:
        fout = os.path.join(settings.PUZZLE_DIR, args.puzzleid + '_' + puzzle.mode + '.out')
    else:
        fout = False
    
    # find solutions
    if args.time:
        [solutions, scores] = optimize_timed(puzzle, args.niter, args.ncool, args.time, draw=args.draw, fout=fout, greedy=args.greedy, start_oligo_conc=args.conc)
    else:
        [solutions, scores] = optimize_n(puzzle, args.niter, args.ncool, args.nsol, draw=args.draw, fout=fout, greedy=args.greedy, start_oligo_conc=args.conc)

if __name__ == '__main__':
    #unittest.main()
    main()


