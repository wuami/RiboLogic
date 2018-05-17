import switch_designer
import sys, os
import json, re
import random
import argparse
import requests
import settings
import varna
import copy
import signal
from collections import OrderedDict

def insert_in_string(str, substr, i):
    return str[0:i] + substr + str[i+len(substr):]

def read_design_from_file(filename, **kwargs):
    """
    reads design information from text file
    """
    inputs = {}
    targets = []
    beginseq = None
    seq_locks = None
    constraints = None
    substr = []
    variables = {}
    
    with open(filename) as f:
        line = f.readline()
        while line:
            # read inputs
            if line.startswith('<'):
                input = f.readline().strip()
                try:
                    inputs[line.strip('<\n')] = {'type':'ligand', 'kD': float(input), 'fold_constraint': f.readline().strip()}
                except:
                    inputs[line.strip('<\n')] = {'type':'RNA', 'sequence':input}
            # read sequence constraints
            elif line.startswith('-sequence'):
                seq = f.readline().strip()
                beginseq = seq.replace('N', 'A')
                seq_locks = ''.join(['o' if c == 'N' else 'x' for c in seq])
                free_positions = [i for i,x in enumerate(seq_locks) if x == 'o']
            # read variable element
            elif line.startswith('-variable'):
                if not beginseq:
                    print 'Must specify sequence before variable element'
                    sys.exit()
                seq = f.readline().strip()
                locations = [m.start() for m in re.finditer('(?=%s)' % 'o'*len(seq), seq_locks)]
                variables[seq] = locations
            # read objectives
            elif line.startswith('>'):
                if not beginseq:
                    print 'Must specify sequence before objectives'
                    sys.exit()
                target = {}
                target['type'] = line.strip('>\n')
                if target['type'] != 'single':
                    target['inputs'] = OrderedDict()
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
                free_positions = [i for i,x in enumerate(target['constrained']) if x == 'o' and i in free_positions]
                if variables:
                    target['variables'] = {}
                    for variable in variables:
                        target['variables'][variable] = {'secstruct': f.readline().strip(),
                                                         'constrained': f.readline().strip()}
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
            elif line.startswith('x'):
                substr.append(line.strip('x\n'))
            line = f.readline()

    free_positions = set(free_positions)
    for seq, pos in variables.items():
        l = len(seq)
        positions = [p for p in pos if set(range(p,p+l)).issubset(free_positions)]
        if 'rpos' in kwargs:
            r = kwargs['rpos']
        else:
            r = random.choice(positions)
        beginseq = insert_in_string(beginseq, seq, r)
        seq_locks = insert_in_string(seq_locks, 'x'*l, r)
        for target in targets:
            target['secstruct'] = insert_in_string(target['secstruct'], target['variables'][seq]['secstruct'], r)
            target['constrained'] = insert_in_string(target['constrained'], target['variables'][seq]['constrained'], r)
    return switch_designer.Design(beginseq, seq_locks, targets, inputs, substrings=substr)

def optimize_n(design, niter, ncool, n, **kwargs):
    # run design n times
    solutions = []
    scores = []
    i = 0 
    attempts = 0
    while i < n:
        design.reset_sequence()
        passkwargs = {key:kwargs[key] for key in ['greedy', 'start_oligo_conc', 'continue_']}
        nfin = design.optimize_sequence(niter, ncool, **passkwargs)
        if design.check_current():
            sol = design.get_solution()
            if sol.sequence not in solutions:
                solutions.append(sol.sequence)
                scores.append(sol.design_score)
                if 'fout' in kwargs and kwargs['fout']:
                    params = ''
                    if 'greedy' in kwargs and kwargs['greedy']:
                        params += 'greedy '
                    with open(kwargs['fout'], 'a') as f:
                        f.write('# %s out of %s iterations, %s coolings, %s\n' % (nfin, niter, ncool, params))
                        f.write('%s\t%1.6f\n' % (sol.sequence, sol.design_score))
                i += 1
                attempts = 0
        else:
            #niter += 500
            print 'best distance: %s' % design.best_design.bp_distance
            print 'final conc: %s' % design.oligo_conc
            attempts += 1
            if attempts == 10:
                break
        print '%s sequence(s) calculated' % i
    return [solutions, scores]

def optimize_timed(design, niter, ncool, time, **kwargs):
    def handler(signum, frame):
        raise Exception('%d elapsed' % time)

    # run design n times
    solutions = []
    scores = []
    niters = []
    i = 0 
    signal.signal(signal.SIGALRM, handler)
    signal.alarm(time)
    try:
        while True:
            design.reset_sequence()
            passkwargs = {key:kwargs[key] for key in ['greedy', 'start_oligo_conc', 'continue_']}
            n = design.optimize_sequence(niter, ncool, **passkwargs)
            if design.check_current():
                sol = design.get_solution()
                solutions.append(sol[0])
                scores.append(sol[2])
                niters.append(n)
                i += 1
            else:
                print 'best distance: %s' % design.best_design.bp_distance
                print 'final conc: %s' % design.oligo_conc
            print '%s sequence(s) calculated' % i
    except Exception, exc:
        print exc
        print '%s sequence(s) calculated in %d seconds' % (i, time)
        print 'average iterations: %d' % (float(sum(niters))/len(niters))
        for i in range(len(solutions)):
            print '\t%s %d %d' % (solutions[i], scores[i], niters[i])
    return [solutions, scores]

def main():
    # parse arguments
    p = argparse.ArgumentParser()
    p.add_argument('filename', help='name of design filename', type=str)
    p.add_argument('-n', '--nsol', help='number of solutions', type=int, default=1)
    p.add_argument('-t', '--time', help='maximum time allowed', type=int)
    p.add_argument('-i', '--niter', help='number of iterations', type=int, default=10000)
    p.add_argument('-o', '--ncool', help='number of times to cool', type=int, default=50)
    p.add_argument('-m', '--mode', help='mode for multi inputs', type=str, default='nupack')
    p.add_argument('-c', '--conc', help='starting oligo concentration', type=float, default=1)
    p.add_argument('--nowrite', help='suppress write to file', default=False, action='store_true')
    p.add_argument('--cont', help='continue optimization after a solution is reached', default=False, action='store_true')
    p.add_argument('--print_', help='print sequences throughout optimization', default=False, action='store_true')
    p.add_argument('--greedy', help='greedy search', default=False, action='store_true')
    p.add_argument('--add_rcs', help='introduce reverse complement of input oligos', default=False, action='store_true')
    args = p.parse_args()

    # read design
    design = read_design_from_file(args.filename)
    if design.default_mode:
        args.mode = design.default_mode
    designer =  switch_designer.SwitchDesigner(os.path.basename(args.filename).split('.')[0], design, **vars(args))
    if not args.nowrite:
        fout = os.path.join(os.path.splitext(args.filename)[0] + '_' + designer.mode + '.out')
    else:
        fout = False
    
    # find solutions
    if args.time:
        [solutions, scores] = optimize_timed(designer, args.niter, args.ncool, args.time, fout=fout, greedy=args.greedy, start_oligo_conc=args.conc, continue_=args.cont)
    else:
        [solutions, scores] = optimize_n(designer, args.niter, args.ncool, args.nsol, fout=fout, greedy=args.greedy, start_oligo_conc=args.conc, continue_=args.cont)

if __name__ == '__main__':
    #unittest.main()
    main()


