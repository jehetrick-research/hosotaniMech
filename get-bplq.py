#!/bin/env python

import sys

fp = open(sys.argv[1], "r")

#print("reading from {}".format(sys.argv[1]))

bplq = {}
nbplq = {}

lines = fp.readlines()
#print(lines)
for L in lines:
    w = L.split()
#    print(w)
    if(len(w)):
        if( w[0] == 'MEASU1'):
            print(w[1], w[len(w)-1])

#            bplq[w[1]] = float(bplq[w[1]]) + float(w[1])
#            nplq[w[1]] += 1

