from functools import wraps
import sys

def intTruth(b):
    if b is True:
        return 1
    return 0

def count(f,lst):
    """Count the number of entries of lst where f(entry) is True"""
    truthList = [intTruth(f(entry)) for entry in lst]
    return sum(truthList)

def hammingDist(a,b):
    #
    d = 0
    for i in xrange(len(a)):
        if a[i] != b[i]:
            d += 1

    return d

    
def deprecated(fn):
    @wraps(fn)
    def wrapper(*args,**kwargs):
        print >>sys.stderr, "Method [%s] is deprecated"%fn.__name__
        return fn(*args,**kwargs)
    return wrapper
    
