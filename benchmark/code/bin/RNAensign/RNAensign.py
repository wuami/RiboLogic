import sys, os
from stochasticSearch import output_conforming_samples

def designWithRNAensign(structure, seed_sequence):
    print len(structure)
    print len(seed_sequence)
    
    samp,nonconsamp = output_conforming_samples(structure,
                                                seed_sequence)
    if samp != []:
        return samp[0]
        

def main(argv):
    structure = argv[0]
    seed_sequence = argv[1]
    ans = designWithRNAensign(structure, seed_sequence)
    if ans is not None:
        print "Sequence: ", ans[0]
        print "Distance: ", ans[1]
        print "Probability: ", ans[2]
    else:
        print "Could not find appropriate sequence" 

if __name__ == '__main__':
    main(sys.argv[1:])

