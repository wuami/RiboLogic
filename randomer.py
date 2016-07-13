import design_utils, fold_utils
import multiprocessing
import string
import numpy as np
import settings



nucleotides = ['A','C','G','U']
def randomer(size):
    return ''.join(np.random.choice(nucleotides) for _ in range(size))

n_targets = 1
sequence = randomer(20)
p = multiprocessing.Pool(n_targets)
result = [None] * n_targets

fold_list = fold_utils.vienna_fold(sequence, bpp=True)
