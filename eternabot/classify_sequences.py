import argparse
import inv_utils
import eterna_utils
import ensemble_design
import design_sequence
import draw_secstructs
import settings, os
import numpy as np
import scipy
from sklearn import cross_validation, naive_bayes, metrics, ensemble, svm, cluster
import pickle

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
        fold_sequence = fold_sequence.replace("&", "")
        secstruct = inv_utils.fold(fold_sequence)[0]
        try:
            design = eterna_utils.get_design_from_sequence(fold_sequence, secstruct)
            bpp_MS2close = [x[2] for x in design['dotplot'] if x[0] == 39 and x[1] == 57]
            if len(bpp_MS2close) == 0:
                bpp_MS2close = [0]
        except:
            design = {'gc':0, 'gu':0, 'ua':0, 'fe':0}
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

def get_top_features(means, stdevs):
    """ get indices of features in descending order of kl divergence """
    kl = []
    for i in range(means.shape[1]):
        m1, m2 = means[:,i]
        s1, s2 = stdevs[:,i]
        kl.append(np.log(s2/s1) + (s1**2+(m1-m2)**2)/(2*s2**2))
    return np.argsort(kl)

def supervised(data, labels, model):
    models = {"SVR": svm.SVR(),
              "RFR": ensemble.RandomForestRegressor(),
              "NB": naive_bayes.GaussianNB()}
    # train and test via 10-fold cross-validation
    kfold = cross_validation.StratifiedKFold(labels, n_folds=10)
    classifier = models[model]
    mean_tpr = 0.0
    mean_fpr = np.linspace(0, 1, 100)
    for i, (train, test) in enumerate(kfold):
        fit = classifier.fit(data[train], labels[train])
        predict = fit.predict_proba(data[test])
        fpr, tpr, thresholds = metrics.roc_curve(labels[test], predict[:,1])
        mean_tpr += scipy.interp(mean_fpr, fpr, tpr)
        mean_tpr[0] = 0.0

    # calculate mean tpr and fpr to generaute AUROC
    mean_tpr /= len(kfold)
    mean_tpr[-1] = 1.0
    mean_auc = metrics.auc(mean_fpr, mean_tpr)
    print mean_auc
    np.savetxt("roc_nb.txt", np.transpose(np.vstack((mean_fpr, mean_tpr))))
    pickle.dump(classifier, open("nb.p",'w'))
    print get_top_features(classifier.theta_, classifier.sigma_)

def unsupervised(data,  method, n_clusters):
    methods = {"KM": cluster.KMeans(n_clusters=n_clusters)}
    
    # cluster using specified method
    clust = methods[method]
    clust.fit(data)
    return clust.labels_

def write_clusters(puzzleid, labels, n_sequences, n_clusters, n_targets):
    # write clusters to html file 
    order = np.argsort(labels)
    breaks = [0]*n_clusters
    for i in labels:
        breaks = [breaks[j]+1 if j >= i else breaks[j] for j in range(n_clusters)]
    draw_secstructs.write_html(puzzleid, n_sequences, n_targets, order, breaks)


def main():
    # parse arugments
    p = argparse.ArgumentParser()
    #p.add_argument("pos", help="filename of positive sequences", type=str)
    #p.add_argument("neg", help="filename of positive sequences", type=str)
    p.add_argument("puzzle", help="filename of puzzle", type=str)
    args = p.parse_args()

    # read in positive and negative training examples
    picklefile = os.path.join(settings.PUZZLE_DIR, args.puzzle + '.p')
    if os.path.isfile(picklefile):
        data, labels = pickle.load(open(picklefile, 'r'))
    else:
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
            puzzle = design_sequence.read_puzzle_json(f.read())

        # create matrices containing features for pos and neg training examples
        data = get_full_matrix_from_sequences(pos_sequences, neg_sequences, puzzle)
        labels = np.array([1]*len(pos_sequences) + [0]*len(neg_sequences))
        pickle.dump([data, labels], open(picklefile, 'w'))

    # perform supervised learning
    #supervised(data, labels, "NB")

    # perform unsupervised learning
    n_clusters = 5
    clusters = unsupervised(data[labels.astype(bool)], "KM", n_clusters)
    write_clusters(args.puzzle, clusters, data.shape[0], n_clusters, 4)

if __name__ == "__main__":
    main()
