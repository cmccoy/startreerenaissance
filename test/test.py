#!/usr/bin/env python

import collections
import itertools
import random

from sys import stderr

from numpy.random import multinomial, seed

seed(0)

# Sequence A: 100 of each base
seqa = ''.join(i * 100 for i in 'ACGT')

print '>A\n{0}'.format(seqa)

print '>B'

B = 0.94
S = 0.02
probs = [[S if i != j else B for j in range(4)]
         for i in range(4)]

counts = ''.join(b * c for p in probs
                 for b, c in zip('ACGT', multinomial(100, p)))
print counts
assert(len(counts) == len(seqa)), (len(seqa), len(counts))
c = collections.Counter(zip(seqa, counts))
for r in 'ACGT':
    print >> stderr, ', '.join(str(c[r, q]) for q in 'ACGT') + ','
