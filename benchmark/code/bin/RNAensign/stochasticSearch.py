import sys,os
import random
from interface import *
from math import exp,log
from utils import *
from config import path_to_rna_fold
from manageTransientFiles import findOutputFile, output_dir

TEMPERATURE = 310.15;
GAS_CONSTANT = 1.985877534 #(Cal K-1 mol-1);
__RT__ = (TEMPERATURE * GAS_CONSTANT)/1000.0;


path_to_rna_fold = '.'

                    
def findAllMFEStructures(listOfSampleSeqs,dangleStr='-d1'):
    samplesFileName = findOutputFile()
    prefix = '> /dev/null'
    seqListToFile(listOfSampleSeqs,samplesFileName,'w',prefix)
    samplesOutputFileName = findOutputFile()
    commandline = '%s/RNAfold -noPS %s -p0 < %s > %s'%(path_to_rna_fold,
                                                       dangleStr,
                                                       samplesFileName,
                                                       samplesOutputFileName)
    os.system(commandline)
    samplesString = open(samplesOutputFileName,'r').read()
    os.remove(samplesFileName)
    os.remove(samplesOutputFileName)
    return mfeprob.findall(samplesString)


def findEntropy(sequence,verbose=False):
    
    d = {}
      
    samplesFileName = findOutputFile()
    name = 'rna' + str(random.randrange(1e8))
    prefix = '> %s'%name
    f = open(samplesFileName,'w')
    f.write(prefix+'\n')
    f.write(sequence+'\n')
    f.close()
    samplesOutputFileName = findOutputFile()
    commandline = '%s/RNAfold -d1 -p < %s > %s'%(path_to_rna_fold,
                                                  samplesFileName,
                                                  samplesOutputFileName)
    os.system(commandline)
    samplesString = open(samplesOutputFileName,'r').read()
    os.remove(samplesFileName)
    os.remove(samplesOutputFileName)
    
    f = open(name + '_dp.ps')
    all_file = f.read()
    second_part = all_file.split('%data starts here\n')[1].split('\n')
    for line in second_part:
        if line.strip() == 'showpage' or 'lbox' in line:
            break
        data = line.split()
        d[(int(data[0]),int(data[1]))] = float(data[2])**2
    entropy = -fsum(d[bp] * log(d[bp],2) for bp in d) / len(sequence)  
    
    if verbose: 
        print '+'.join([str(d.get(bp,1e-6) * log(d.get(bp,1e-6),2)) for bp in gen_basePairs(structure)])

    os.remove( name + '_dp.ps')
    os.remove( name + '_ss.ps')
        
    
    assert entropy >=0 
    return entropy


def mfeStructProb(sequence):
    return findAllMFEStructures([sequence])[0][2]


def output_conforming_samples(targetStruct,sequence,sampleNumber = 1000,
                            maxNumMutations = None): 
    """
    Run RNAensign, and output samples that achieve a given structure
    as MFE structure.

    Input:
    <targetStruct>
    <sequence>: the seed sequence
    <sampleNumber>: Number of samples for each mutation number
    <maxNumMutations> (default None): the maximum number of mutations 
    that RNAensign will perform.  If this is None, it can perform any
    number of mutations.  By decreasing this, we can substantially
    decrease the running time, but may not find conforming sequences if
    more mutations are needed in order to do so.

    Output:
    The list of samples (or [] if there were none).  Each sample is a
    tuple consisting of sequence, distance from seed, and probability of
    achieving the target structure.
    """
    sampleslist = []
    nonconformingsampleslist = []
    moreopts = ''
    paramFileName = findOutputFile()
    f = open(paramFileName,'w')
    if maxNumMutations is None or maxNumMutations > len(sequence):
        maxNumMutations = len(sequence)
    for i in xrange(maxNumMutations):
        f.write('%s %s\n'%(i+1,sampleNumber))
    
    f.close()
    
    samples,superoptIter,\
        energiesDict = generateParametrizedSamples(targetStruct, sequence, 
                                                   paramFileName,
                                                   True,True,moreopts)
    os.remove(paramFileName)
    for samplesWithKMut in samples:
        for sampseq, mfestruct, p in findAllMFEStructures(samplesWithKMut):
            assert len(sampseq) == len(sequence)
            if len(sampseq) != len(sequence): 
                continue
            if mfestruct == targetStruct:
                sampleslist.append((sampseq,int(hammingDist(sampseq,sequence))
                                    ,float(p)))
    
            else:
                nonconformingsampleslist.append((sampseq,int(hammingDist(sampseq,sequence)),float(p)))
                
    try:
        os.rmdir(output_dir)
    except OSError:
        pass

    return sampleslist,nonconformingsampleslist

