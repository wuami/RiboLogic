import os
import pandas as pd
import numpy as np
import argparse
from Levenshtein import distance
import multiprocessing as mp
from functools import partial


def get_pairmap(secstruct):
  """Get the pairing partners of each base in a secondary structure.

  List of indices of pairing partners for each base, with -1 if unpaired.

  Args:
    secstruct: secondary structure in dot-bracket notation

  Returns:
    list of pairing partners
  """
  pairs = [-1] * len(secstruct)
  bases = []
  for i, base in enumerate(secstruct):
    if base == '(':
      bases.append(i)
    elif base == ')':
      j = bases.pop()
      pairs[i] = j
      pairs[j] = i
  return pairs

def bp_distance(secstruct1, secstruct2):
  """Get base pair distance between two secondary structures.

  This computes the number of base pairs that need to be broken or formed
  to change one secondary structure to the other.

  Args:
    secstruct1: first secondary structure in dot bracket notation
    secstruct2: second secondary structure in dot bracket notation

  Returns:
    int: base pair distance value
  """
  assert len(secstruct1) == len(secstruct2), (secstruct1, secstruct2)
  pm1 = get_pairmap(secstruct1)
  pm2 = get_pairmap(secstruct2)
  dist = 0
  for i in range(len(pm1)):
    if pm1[i] != pm2[i] and (pm1[i] > i or pm2[i] > i):
      dist += 1
  return dist


def get_distances(values, function=distance):
  """Get Levenshtein distance between last value and all previous.

  Args:
    values: list of strings

  Returns:
    list of distances between values 0..n-2 and value n-1
  """
  n = len(values)
  return np.array([function(values[i], values[n - 1]) for i in range(n - 1)])


def main():
  p = argparse.ArgumentParser()
  p.add_argument('filename', help='name of input file')
  p.add_argument('column', help='name of column in input file to compute '
                 'distances on')
  args = p.parse_args()

  values = pd.read_csv(args.filename, sep='\t')[args.column]
  n = len(values)

  p = mp.Pool()
  # distance metric is determined based on column name
  result = map(partial(get_distances, function=bp_distance if 'secstruct' in args.column else distance), [values[0:i+1] for i in range(n)])
  result = np.array([np.pad(row, (0, n - len(row)), 'constant')
                     for row in result])
  np.savetxt('%s.%s.dist' % (os.path.splitext(args.filename)[0], args.column), result, '%d',
             delimiter='\t')


if __name__ == '__main__':
  main()
