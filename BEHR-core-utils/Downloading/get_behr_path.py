from __future__ import print_function

import os
import re
import sys

import pdb

_util_root = os.path.abspath(os.path.realpath(os.path.join(os.path.dirname(__file__), '..')))
_paths_file = os.path.join(_util_root, 'Utils', 'Constants', 'behr_paths.m')

def get_path(pathname):
    regex = re.compile('{}.*='.format(pathname))
    with open(_paths_file, 'r') as fobj:
        for line in fobj:
            if regex.search(line) is not None:
                m = re.search("(?<=').+(?=')", line)
                return m.group()

def main():
    pathname = sys.argv[1]
    print(get_path(pathname))
    exit(0)

if __name__ == '__main__':
    main()
