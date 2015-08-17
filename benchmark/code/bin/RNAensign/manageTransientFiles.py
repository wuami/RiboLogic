import sys,os,pickle,tempfile
import random
from time import strftime

constraintDictFileName = 'constraintsDict.txt'

output_dir = tempfile.mkdtemp(prefix='nn-tmp-')

def timeStamp( includeRand = None ):
    stamp = strftime('%Y%m%d-%H%M%S')
    if includeRand is not None:
        stamp += '_'+str(random.randint(0,1e6))
        
    return stamp
        
def handleStructure(structure):
    constraintDictFile = open(constraintDictFileName,'r')
    cfDict = pickle.load(constraintDictFile)
    #dictionary mapping structures to constraint files
    constraintDictFile.close()
    
    try:
        constraintFileName = cfDict[structure]
    except KeyError:
        constraintDictFile = open(constraintDictFileName,'w')
        constraintFileName = strftime('constraints%Y%m%d-%H%m%s.cst')
        f = open(constraintFileName,'w')
        f.write(structure + '\n')
        f.close()
        cfDict[structure] = constraintFileName
        pickle.dump(constraintFileDict,constraintDictFile)
        constraintDictFile.close()
    return constraintFileName
        
def structToAlpha( struct ):
    return struct.replace('(','L').replace(')','R').replace('.','D')

def findConstraintFile(struct):
    alphaStruct = structToAlpha(struct)
    structFileName = 'constraints/constraint_' + alphaStruct + '.cst'
    try:
        f = open(structFileName,'r')
        structInFile = f.readline()
        if structInFile.rstrip('\n') != struct:
            raise IOError
    
    except IOError:
        f = open(structFileName,'w')
        f.write(struct+'\n')

    return structFileName

def findOutputFile():
    try:
        subdir = os.environ['HOSTNAME']
    except KeyError:
        subdir = ''
    if 'EXPERIMENT' in os.environ: subdir += '_' + os.environ['EXPERIMENT']
    try: os.makedirs('output/%s' % subdir)
    except OSError, e: pass # already exists
    return '%s/output_%s.txt'%(output_dir,timeStamp( 1 ))


#Necessary for NUPACK
#all the files it uses are of the form PREFIX.init, PREFIX.out, etc
#for the same prefix
def makeNUPACKFilePrefix():
#TODO This is copied code
    nupack_output_dir = tempfile.mkdtemp(prefix='nupack-tmp')#'/tmp/nupack-tmp'

    try:
        subdir = os.environ['HOSTNAME']
    except KeyError:
        subdir = ''
    if 'EXPERIMENT' in os.environ: subdir += '_' + os.environ['EXPERIMENT']
    try: os.makedirs('output/%s' % subdir)
    except OSError, e: pass # already exists
    return '%s/prefix_%s'%(nupack_output_dir,timeStamp( 1 )),nupack_output_dir
    
