import argparse
import inv_utils
import eterna_utils
import ensemble_design
import switch_design
import settings, os
import numpy as np
from sklearn import cross_validation, naive_bayes, metrics

def get_random_sequences(n, length):
    """ generate n random sequences of length specified """
    bases = "AUGC"
    sequences = []
    for i in range(n):
        r = np.random.randint(4, size=length)
        sequence = [bases[i] for i in r]
        sequences.append(ensemble_design.get_sequence_string(sequence))
    return sequences

def get_vector_from_sequence(sequence, puzzle):
    """ create vector to represent a sequence """
    vector = []

    for target in puzzle.targets:
        fold_sequence = puzzle.get_fold_sequence(sequence, target)
        fold_sequence = fold_sequence.replace("&", "\&")
        secstruct = inv_utils.fold(fold_sequence)
        design = eterna_utils.get_design_from_sequence(fold_sequence, secstruct)
        bpp_MS2close = [x[2] for x in design['dotplot'] if x[0] == 39 and x[1] == 57]
        if len(bpp_MS2close) == 0:
            bpp_MS2close = [0]
        vector.extend([design['gc'], design['gu'], design['ua'], design['fe']])#, bpp_MS2close[0]])
    
    a = t = c = g = 0
    for i in range(len(sequence)):
        if sequence[i] == "A":
            a += 1
        elif sequence[i] == "T":
            t += 1
        elif sequence[i] == "C":
            c += 1
        elif sequence[i] == "G":
            g += 1
    vector.extend([a, t, c, g])

    return vector

def get_matrix_from_sequences(sequences, puzzle):
    """ create a matrix to represent a list of sequences """
    matrix = []
    for sequence in sequences:
        matrix.append(get_vector_from_sequence(sequence, puzzle))
    return np.array(matrix)

def get_full_matrix_from_sequences(pos_sequences, neg_sequences, puzzle):
    """ create combined matrix for pos and neg sequences """
    return np.vstack((get_matrix_from_sequences(pos_sequences, puzzle), get_matrix_from_sequences(neg_sequences, puzzle)))

def main():
    # parse arugments
    p = argparse.ArgumentParser()
    #p.add_argument("pos", help="filename of positive sequences", type=str)
    #p.add_argument("neg", help="filename of positive sequences", type=str)
    p.add_argument("puzzle", help="filename of puzzle", type=str)
    args = p.parse_args()

    # read in positive and negative training examples
    pos_sequences = []
    with open(os.path.join(settings.PUZZLE_DIR, args.puzzle + '.out'), 'r') as f:
        for line in f:
            if not line.startswith('#'):
                pos_sequences.append(line.split()[0])
    neg_sequences = get_random_sequences(len(pos_sequences), len(pos_sequences[0]))
    #neg_sequences = []
    #with open(args.neg, 'r') as f:
    #    for line in f:
    #        if not line.beginswith('#'):
    #            neg_sequences.append(line.split()[0])

    # get puzzle
    with open(os.path.join(settings.PUZZLE_DIR, args.puzzle + '.json'), 'r') as f:
        puzzle = switch_design.read_puzzle_json(f.read())

    # create matrices containing features for pos and neg training examples
    data = get_full_matrix_from_sequences(pos_sequences, neg_sequences, puzzle)
    labels = np.array([1]*len(pos_sequences) + [0]*len(neg_sequences))

    # train and test
    data_train, data_test, labels_train, labels_test = cross_validation.train_test_split(data, labels, test_size=0.3)
    clf = naive_bayes.GaussianNB()
    clf.fit(data_train, labels_train)
    print labels_test, clf.predict_proba(data_test)
    fpr, tpr, thresholds = metrics.roc_curve(labels_test, clf.predict_proba(data_test)[:,1])
    print metrics.auc(fpr, tpr)
    np.savetxt("roc_nb.txt", np.vstack((fpr, tpr)))

if __name__ == "__main__":
    main()
