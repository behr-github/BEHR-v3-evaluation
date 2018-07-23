#!/usr/bin/python
__author__ = 'Josh'
try:
    import cPickle as pickle
except ImportError:
    import pickle

def make_pickle(obj, filename):
    with open(filename, 'wb') as f:
        pickle.dump(obj, f, pickle.HIGHEST_PROTOCOL)
