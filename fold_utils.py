import subprocess, os, sys, settings
import re, string, random, itertools
import numpy as np
import warnings

paramfile = 'rna'

def vienna_fold(sequence, constraint=False, bpp=False, version=settings.vienna_version):
    """
    folds sequence using Vienna

    args:
    seq is the sequence string

    returns:
    secondary structure
    """
    if version in ['1.8.5', '1.8', '1']:
        if settings.os_ == 'osx':
            print 'Vienna version 1.8.5 not supported on OSX'
            sys.exit()
        version = '1'
    elif version in ['2.1.9', '2.1', '2']:
        version = '2'
    else:
        print 'Vienna version must be 1.8.5 or 2.1.9'
        sys.exit()
    filename = '%s/%s' % (settings.TEMP_DIR, ''.join(random.sample(string.lowercase,5)))
    with open('%s.fa' % filename, 'w') as f:
        f.write('>%s\n' % os.path.basename(filename))
        f.write('%s\n' % sequence)
        if constraint:
            options = ' -C'
            f.write('%s\n' % constraint)
        else:
            options = ''
    if bpp:
        options += ' -p'
    else:
        if version == '1':
            options += ' -noPS'
        else:
            options += ' --noPS'
    if '&' in sequence:
        if sequence.count('&') > 1:
            print 'Cannot handle more than 2 strands with Vienna - try the nupack option'
            sys.exit()
        command = 'RNAcofold' + version
    else:
        command = 'RNAfold' + version
    output = subprocess.check_output(os.path.join(settings.VIENNA_DIR,command) + options + ' -T 37.0 < %s' % filename + '.fa', shell=True)

    # parse the result
    toks = re.search('([AUGC]+)\s*([\.\)\(&]+)\s+\(\s*([-0-9\.]+)\s*\)', output)

    if bpp:
        bpp_matrix = np.zeros((len(sequence), len(sequence)))
        # get info from output file
        psfile = '%s_dp.ps' % os.path.basename(filename)
        if os.path.isfile(psfile):
            with open(psfile) as f:
                ps = f.read()

            # create bpp matrix from file
            lines = re.findall('(\d+)\s+(\d+)\s+(\d*\.*\d*)\s+ubox',ps)
            for ii in range(0,len(lines)):
                bpp_matrix[int(lines[ii][0])-1, int(lines[ii][1])-1] = float(lines[ii][2])
        else:
            warnings.warn('dotplot file %s_dp.ps not found' % filename)
        os.system('rm %s*' % filename)
        os.system('rm %s*' % os.path.basename(filename))
        return [toks.group(2), float(toks.group(3)), np.array(bpp_matrix)]
    
    os.system('rm %s*' % filename)
    return [toks.group(2), float(toks.group(3))]

def get_orderings(n):
    """
    get all possible orderings including last strand
    """
    all = []
    # loop over number of strands
    for i in range(1,n):
        # loop over each possible combination
        for order in list(itertools.combinations(range(1,n), i)):
            # add last strand at each possible position
            for j in range(1,i+1):
                order_list = list(order)
                order_list.insert(j,n)
                all.append(order_list)
    return all

def nupack_fold(seq, oligo_conc=1, bpp=False):
    """
    finds most prevalent structure using nupack partition function
    """
    if '&' in seq:
        return nupack_fold_multi(seq, oligo_conc, bpp)
    else:
        return nupack_fold_single(seq, bpp)

def nupack_fold_single(seq, bpp=False):
    """
    finds most prevalent structure using nupack partition function for single strand
    """
    rand_string = ''.join(random.choice(string.ascii_lowercase) for _ in range(6))
    filename = '%s/%s' % (settings.TEMP_DIR, rand_string)
    options = ['-material', paramfile]
    with open('%s.in' % filename, 'w') as f:
        f.write('%s\n' % seq)
    p = subprocess.Popen([os.path.join(settings.NUPACK_DIR,'mfe')] + options + [filename], stdout=subprocess.PIPE, stdin=subprocess.PIPE, stderr=subprocess.PIPE)
    stdout, stderr = p.communicate()
    if p.returncode != 0:
        print stdout
        print stderr
        raise ValueError('mfe command failed for %s' % rand_string)
    result = ['.'*len(seq), 0, [1]]
    with open('%s.mfe' % filename) as f:
        line = f.readline()
        while line:
            if not line.startswith('%') and line.strip():
                energy = float(f.readline().strip())
                secstruct = f.readline().strip()
                result = [secstruct, energy, [1]]
                break
            line = f.readline()
    if bpp:
        p = subprocess.Popen([os.path.join(settings.NUPACK_DIR,'pairs')] + options + [filename], stdout=subprocess.PIPE, stdin=subprocess.PIPE, stderr=subprocess.PIPE)
        stdout, stderr = p.communicate()
        bpp_matrix = []
        if p.returncode != 0:
            print stdout
            print stderr
            raise ValueError('mfe command failed for %s' % rand_string)
        with open('%s.ppairs' % filename) as f:
            line = f.readline()
            while line:
                if not line.startswith('%') and line.strip():
                    bpp_matrix = nupack_read_bpp(f, len(seq))
                line = f.readline()
        result.append(bpp_matrix)
    os.system('rm %s*' % filename)
    return result

def nupack_fold_multi(seq, oligo_conc=1, bpp=False):
    """
    finds most prevalent structure using nupack partition function for multistrand
    """
    rand_string = ''.join(random.choice(string.ascii_lowercase) for _ in range(6))
    filename = '%s/%s' % (settings.TEMP_DIR, rand_string)
    split = seq.split('&')
    with open('%s.in' % filename, 'w') as f:
        f.write('%s\n' % len(split))
        for s in split:
            f.write('%s\n' % s)
        f.write('1\n')
    orderings = get_orderings(len(split))
    with open('%s.list' % filename, 'w') as f:
        for ordering in orderings:
            f.write('%s\n' % ' '.join([str(x) for x in ordering]))
    with open('%s.con' % filename, 'w') as f:
        if isinstance(oligo_conc, list):
            assert len(split) == len(oligo_conc) + 1, \
                'length of concentrations must be one less than number of strands'
            f.write('%s\n' % '\n'.join([str(x) for x in oligo_conc]))
        else:
            f.write('%s\n' % oligo_conc * (len(split)-1))
        f.write('1e-9\n')
    options = ['-material', paramfile, '-ordered', '-mfe']#, '-quiet']
    if bpp:
        options.append('-pairs')
    p = subprocess.Popen([os.path.join(settings.NUPACK_DIR,'complexes')] + options + [filename], stdout=subprocess.PIPE, stdin=subprocess.PIPE, stderr=subprocess.PIPE)
    stdout, stderr = p.communicate()
    if p.returncode != 0:
        print stdout
        print stderr
        raise ValueError('complexes command failed for %s' % rand_string)
    p = subprocess.Popen([os.path.join(settings.NUPACK_DIR,'concentrations'), '-ordered', filename], stdout=subprocess.PIPE, stdin=subprocess.PIPE, stderr=subprocess.PIPE)
    stdout, stderr = p.communicate()
    if p.returncode != 0:
        print stdout
        print stderr
        raise ValueError('concentrations command failed for %s' % rand_string)
    # get mfe
    complex = False
    with open('%s.eq' % filename) as f_eq:
        for line in f_eq:
            if not line.startswith('%'):
                complex = line.strip().split()
                if float(complex[len(split)+1]):
                    break
    if not complex:
        os.system('rm %s*' % filename)
        if bpp:
            return ['.'*len(seq), 0, range(1,len(split)+1), []]
        return ['.'*len(seq), 0, range(1,len(split)+1)]
    # get strand ordering
    with open('%s.ocx-key' % filename) as f_key:
        for line in f_key:
            if line.startswith('%s\t%s' % (complex[0], complex[1])):
                strands = [int(x) for x in line.strip().split()[2:]]
                break
    # get mfe structure
    with open('%s.ocx-mfe' % filename) as f_mfe:
        line = f_mfe.readline()
        while line:
            if line.startswith('%% complex%s-order%s' % (complex[0], complex[1])):
                f_mfe.readline()
                energy = f_mfe.readline().strip()
                secstruct = f_mfe.readline().strip()
                break
            line = f_mfe.readline()

    # get full secondary structure
    for i in range(len(split)):
        if i+1 not in strands:
            secstruct += '&' + '.'*len(split[i])
            strands.append(i+1)
    
    if bpp:
        bpp_matrix = []
        with open('%s.cx-epairs' % filename) as f_pairs:
            line = f_pairs.readline()
            while line:
                if line.startswith('%% complex%s' % complex[0]):
                    bpp_matrix = nupack_read_bpp(f_pairs, len(seq))
                    break
                line = f_pairs.readline()
        os.system('rm %s*' % filename)
        return [secstruct.replace('+', '&'), float(energy), strands, bpp_matrix]

    os.system('rm %s*' % filename)
    return [secstruct.replace('+', '&'), float(energy), strands]

def nupack_read_bpp(f, n):
    """
    read a base pair probability matrix from a nupack output file
    """
    bpp_matrix = np.zeros((n,n))
    f.readline()
    line = f.readline()
    while line and not line.startswith('%'):
        bp = line.strip().split()
        if int(bp[1])-1 != n:
            bpp_matrix[int(bp[0])-1, int(bp[1])-1] = float(bp[2])
        line = f.readline()
    return bpp_matrix


def nupack_energy(seq, secstruct):
    rand_string = ''.join(random.choice(string.ascii_lowercase) for _ in range(6))
    filename = '%s/%s' % (settings.TEMP_DIR, rand_string)
    multi = '&' in seq
    split = seq.split('&')
    with open('%s.in' % filename, 'w') as f:
        if multi:
            f.write('%s\n' % len(split))
        for seq in split:
            f.write('%s\n' % seq)
        if multi:
            f.write('%s\n' % ' '.join([str(i) for i in secstruct[1]]))
            f.write(secstruct[0].replace('&', '+'))
        else:
            f.write(secstruct.replace('&', '+'))
    if '&' in seq:
        options = ['multi']
    else:
        options = []
    result = subprocess.check_output([os.path.join(settings.NUPACK_DIR,'energy')] + options + [filename])
    os.system('rm %s*' % filename)
    return float(result.strip().split('\n')[-1])


