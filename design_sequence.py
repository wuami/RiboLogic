import switch_designer
import os
import re
import argparse
import signal

def optimize_n(design, niter, ncool, n, maxattempts=1, **kwargs):
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
            if attempts == maxattempts:
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
    design = switch_designer.read_design_from_file(args.filename)
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


