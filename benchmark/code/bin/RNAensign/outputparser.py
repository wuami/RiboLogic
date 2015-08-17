import sys, os, re

optStructures = re.compile(r"""
>\s(?P<numMutations>\d+).*\n
(?P<sequence>[a-zA-Z]+)\s+[(](?P<energy>[.\d-]+)[)]\n
(?P<structure>[(.)]+)
    """,re.VERBOSE)

sample = re.compile(r"""
(?P<sequence>[AaUuGgCc]+)\s+[(](?P<energy>[.\d-]+)[)]\n
(?P<structure>[(.)]+)
""",re.VERBOSE)

mfeprob = re.compile(r"""
(?P<sequence>[AaUuGgCc]+)\n
(?P<mfestruct>[(.)]+).+\n
.+\n
.+(?P<prob>\d+[.]\d+)[;]
""",re.VERBOSE)
#parse a run from RNAfold


ensdivers = re.compile(r"""
diversity\s(\d+[.]\d+)
""",re.VERBOSE)

RNASSDwebserverresult = re.compile(r"""
5'-(?P<sequence>[AaUuGgCc]+)-3'
""",re.VERBOSE)

INFORNAwebserverresult = re.compile(r"""
[<]TEXTAREA [^\n]*[>]\n
(?P<sequence>[AaUuGgCc]+)\n
##[<]/TEXTAREA[>]
""",re.VERBOSE)


def parseMutantRun( mutantData ):
    return optStructures.finditer( mutantData )
    
def parseSamples( mutantData):
    l = mutantData.split('>>')
    return sample.finditer(l[-1])
    #in order to only take the correct portion of the file

def parseParametrizedSamples( mutantData ):
    smp = mutantData.split('>>')[-1]
    itOfSamples = [sample.finditer(block) for block in smp.split('>')] 
    return itOfSamples

def getPartitionFunctionValue(mutantData, i = 0):
    linesOfMutantData = mutantData.split('\n')
    return float(linesOfMutantData[i].split()[1])
        
def getSequenceFromRNASSDResult(resultstring):
    return RNASSDwebserverresult.search(resultstring).group('sequence')

def getSequenceFromINFORNAResult(resultstring):
    return INFORNAwebserverresult.search(resultstring).group('sequence')


###################################
#parse NUPACK output
