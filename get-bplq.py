#!/bin/env python

import sys

fp = open(sys.argv[0], "r")


lines = fp.readlines()
for L in lines:
    w = L.split()
    print(w)
#    if(len(w)):
    if( w[0] == 'MEASU1'):
       print(w[1], w[len(w)-1])
 
