#!/bin/env python

# JH 30 Sept 2016
# Generate in.file for ks_spectrum program
#
# measure PsiBar-Psi on lattices in the (beta, h) plane
#
# Hardcoded parameters:
#
#    nx,...,nt = 16^3 x 4
#    steps_per_trajectory = 1
#    no_gauge_fixing
#    lattice name stub "l164bxxxhxx.n"
#
# Input parameters:
#                     default
#    beta_min,max,delta
#    h_min,max,delta
#    output_infile  -> in.pbploop

output_infile = 'in.pbploop'

import sys
def eprint(*args, **kwargs):
    print(*args, file=sys.stderr, **kwargs)

beta_min = float(input("beta_min = "))
beta_max = float(input("beta_max = "))
beta_delta = float(input("beta_delta = "))
h_min = float(input("h_min = "))
h_max = float(input("h_max = "))
h_delta = float(input("h_delta = "))

prompt = "output in.file = [{}]: ".format(output_infile)
t = input(prompt)
output_infile = t if t != '' else output_infile

eprint("beta: {} to {} by {}".format(beta_min, beta_max, beta_delta))
eprint("h   : {} to {} by {}".format(h_min, h_max, h_delta))
eprint("output written to: {}".format(output_infile))

ofile = open(output_infile, "w")

import random
header = """prompt 0
nx 16
ny 16
nz 16
nt 4
iseed {}
job_id 777

""".format(random.randint(0,10000))

print(header, file=ofile)


# Read in u0
import pandas as pd
u0df = pd.read_csv("u0list.txt", sep=" ", names=["startfile", "beta", "h", "ssplaq", 'stplaq', 'u0'])


# loop over lattices
import numpy as np
nb = round((beta_max - beta_min)/beta_delta + 1)
nh = round((h_max - h_min)/h_delta + 1)
barray = np.linspace(beta_min, beta_max, nb)
#barray = barray[::-1]

for h in np.linspace(h_min, h_max, nh):
#    barray = barray[::-1]
    for b in barray:
        startlat = "l164b{:3d}h{:02d}".format(int(round(100*b)),
                                             int(round(100*h)))
        u0 = float(u0df[u0df.startfile==startlat].u0)

        stanza = """reload_serial lattices/{}
u0 {}
no_gauge_fix
forget
staple_weight 0
ape_iter 0
coordinate_origin 0 0 0 0

# Chiral condensate and related measurements

number_of_pbp_masses 2

max_cg_iterations 300
max_cg_restarts 5
npbp_reps 8
prec_pbp 1

mass 0.02
error_for_propagator 1e-6
rel_error_for_propagator 0
mass 0.04
error_for_propagator 1e-6
rel_error_for_propagator 0

number_of_base_sources 0 number_of_modified_sources 0 number_of_sets 0
number_of_quarks 0 number_of_mesons 0 number_of_baryons 0
""".format(startlat, u0)

        print(stanza, file=ofile)





