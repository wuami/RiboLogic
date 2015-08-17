import sys,os,re
import pickle
from time import strftime
import random
from math import *
import subprocess
import tempfile

from outputparser import *
from config import path_to_rna_mutants

_globalJobID = None


# Useful for limiting the number of mutations
# to perform.
def findmaxk(filename):
    f=open(filename,'r');
    linelist=f.readlines();
    f.close();

    maxk=0;
    for line in linelist:
        field=line.split(' ');
        k=int(field[0]);
        if k>maxk:
            maxk=k;

    return maxk;

class RNAState(object):
    def __init__(self,sequence,structure,numMutations):
        self.sequence = sequence.upper()
        self.structure = structure
        self.numMutations = int(numMutations)
        #should we include mutation number?  In that case, it should
        #be typecast to an integer


def _runRNAMutants(sequence,structure = None, extraOpts = '',isSample = False,
                   isParamSample = False, displayOutput = False,
                   returnRunString = False):
    """Generic function for running the RNAmutants software.  
    We generate the lowest-energy mutants starting at [sequence]; 
    [structure] denotes structure constraints (no constraints is the default), 
    and [extraOpts] are extra options.  

    The <isSample> flag, when set to True, tells us to output samples 
    rather than mutation profile.  The <isParamSample> flag
    tells us to output a list of iterators, one for each
    line in the parameter file.  For example, if the parameter file
    is:
    --------------------
       1 1000
       2 500
       3 1000
    --------------------
    then the first entry of the list is an iterator for 1000 samples
    of 1-mutants, the second entry is an iterator for 500 samples of 2-mutants,
    and the third entry is for 

    Note that we will have to specify the parameter file (and the 
    extra options that denote parametrized sampling) in functions 
    calling _runRNAMutants (see generateParametrizedSamples, and 
    the setting of extraOpts therein). 
    """    

    # setup directories
    inputdir = tempfile.mkdtemp(prefix="nn-rm-input-");
    outputdir = tempfile.mkdtemp(prefix="nn-rm-output-");
    assert os.path.exists(inputdir)
    assert os.path.exists(outputdir)

    # setup filenames
    validfiles=False;
    while not validfiles:
        myhashcode = "%d" % (random.randint(0,1e9));
        inputfile  = "%s/input_%s.txt"  % (inputdir,myhashcode);
        outputfile = "%s/output_%s.txt" % (outputdir,myhashcode);
        if not os.path.isfile(inputfile) and not os.path.isfile(outputfile):
            validfiles = True;

    # touch output filename
    f = open(outputfile,'w')
    f.write('')
    f.close()

    # extract path to RNAmutants	 
    path = path_to_rna_mutants
    #remove last front slash in path
    path = path.strip().rstrip('/')

    # create command line 
    commandLine = '%s/RNAmutants --library %s/lib --file %s ' % (path,path,inputfile);

    ff = open(inputfile,'w');
    ff.write(sequence + '\n');
    if structure is not None:
        ff.write(sequence.upper().replace('A','o').replace('C','o').replace('G','o').replace('U','o').replace('T','o') + '\n');
        ff.write(structure + '\n');
        commandLine += ' -C ';			
    ff.close();
   
    if getJobID() is not None:
        commandLine += "--job-id '%s' " % getJobID()
    commandLine += extraOpts + ' '
        
    commandLine+='> %s'%outputfile 
    try:
        subprocess.call(commandLine, shell=True)
    except KeyboardInterrupt:
        os.remove(inputfile)
        os.remove(outputfile)
        os.rmdir(inputdir)
        os.rmdir(outputdir)

        raise
        
        #sys.exit()
        
    #just return the output on the screen without printing it to a file
    if displayOutput:
        print commandLine
        return None
    
    # read RNAmutants output
    mutantData = open(outputfile,'r').read()
    
   # delete the temporary files in the output 
    os.remove(outputfile);
    os.remove(inputfile)
    #delete temporary directories
    os.rmdir(inputdir)
    os.rmdir(outputdir)

    ## extract data and return iterators
    
    if returnRunString:
        return mutantData


    iterMFE = parseMutantRun(mutantData)

    #use different functions from outputparser.py depending on [isSample] flag
    iterSMP = None #default
    if isSample:
        iterSMP = parseSamples(mutantData)
    elif isParamSample:
        iterSMP = parseParametrizedSamples( mutantData )
        # careful: isSample and isParamSample should not be simultaneously true
    
    iter = parseMutantRun(mutantData)
            
    return iter, iterSMP;


def findMFEStructure(sequence):
    iter, iterSMP = _runRNAMutants(sequence,extraOpts  = '--mutations 0')
    return iter.next().group('structure')


def generateSamples(structure,seedSequence,numSamples,extraOpts = ''):
    extraOpts += '-n %s'%numSamples
    iterSuperopt, iter = _runRNAMutants(seedSequence, structure, 
                                               extraOpts,
                                               isSample = True)
    samples = []
    for mutant in iter:
        samples.append(mutant.group('sequence').upper())
    return samples
    
    
def generateParametrizedSamples(structure, seedSequence, paramFileName,
                                returnSuperOptIter = False,
                                returnEnergies = False,
                                moreoptions=''):
    
    """Generate samples consistent
    with <structure> and whose mutation profile
    from <seedSequence> is parametrized by <paramFileName>""" 
    
    maxk = findmaxk(paramFileName);
    #Necessary so that in the call to RNAmutants, you 
    #only do up to <maxk> mutations, thus saving time

    extraOpts = '--sample-file %s --mutations %d '% (paramFileName,maxk);
    extraOpts += moreoptions
    try:
        iterSuperopt, listOfIters = _runRNAMutants(seedSequence, structure, 
                                               extraOpts,
                                               isParamSample = True)
    except KeyboardInterrupt:
        os.remove(paramFileName)
        dir = paramFileName.rsplit('/',1)[0]
        print >>sys.stderr, dir
        os.rmdir(dir)
        sys.exit()
    #a list of iterators, one for each line in the parameter file
    
    #now we turn it into a list of sample sequences; each 
    #entry in the list will itself be a list
    #of sequences, corresponding to one line in the parameter file.

    listOfSamples = []
    energyDict = {}
    for iter in listOfIters:
        samplesWithKMutants = []
        for mutant in iter:
            samp_sequence = mutant.group('sequence').upper()
            energyDict[samp_sequence] = float(mutant.group('energy'))
            samplesWithKMutants.append(samp_sequence)
        listOfSamples.append(samplesWithKMutants)
    if returnSuperOptIter:
        if returnEnergies:
            return listOfSamples[1:],iterSuperopt,energyDict
        return listOfSamples[1:],iterSuperopt
    return listOfSamples[1:]
	#for some reason the first element of listOfSamples is the empty list
	#so we use this slicing

        

        
def sequenceStatistics(listOfSequences):
    """Gather statistics on nucleotide frequencies in the 
    samples <listOfSequences>.  We think of <listOfSequences>
    as drawn from from high temperature ensemble 
    consistent with a given structure.  
    Hope to eventually feed this to GenRGenS."""

    stats = {}
    fullString = "".join(listOfSequences)
    totalLength = len(fullString)
    for nucleotide in ['A','U','G','C']:
        stats[nucleotide] = float(fullString.count(nucleotide))/totalLength
   
    return stats


def PurineProportion(seq):
    if len(seq) != 0:
        return float( seq.count('A') + seq.count('G')) / len(seq)
    return 'N/A'

def GCProportion(seq):
    if len(seq) != 0:
        return float( seq.count('G') + seq.count('C')) / len(seq)
    return 'N/A'
    
def PurineProportionOfList(listOfSeqs):
    return GCProportion("".join(listOfSeqs))

def GCProportionOfList(listOfSeqs):
    return GCProportion("".join(listOfSeqs))

def fileToSeqList(fileName):
    seqList = []
    file = open(fileName,'r')
    for line in file:
        seqList.append(line.replace(" ","").replace('\n',''))
    file.close()
    return seqList	 	

def seqListToFile(seqList,fileName,method = 'a',prefix = None):
    file = open(fileName,method)
    if prefix is not None:
        file.write(prefix+'\n')
    for seq in seqList:
        file.write( seq +'\n')
    file.close()

def setJobID(id):
    global _globalJobID
    if type(id) is not str: id = '%s' % id
    #_globalJobID = id if len(id) > 0 else None
    # incompatible with earlier Python versions
    if len(id) > 0: 
        _globalJobID = id
    else: 
        _globalJobID = None

def getJobID():
    return _globalJobID

