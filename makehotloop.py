#!/bin/env python

# JH 30 Sept 2016
# Generate in.file for su3_trP program
#
# Use to build an array of lattices in the (beta, h) plane
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
#    initial warms  -> 50000
#    warms btw lats ->  5000
#    trajecs        ->  1000 (i.e. traj_btw_meas=10 -> 100 measurements)
#    traj_between_meas -> 10
#    initial_lat    ->  l164b560.20000
#    output_infile  -> in.loop

import sys
def eprint(*args, **kwargs):
    print(*args, file=sys.stderr, **kwargs)


initial_warms = 50000
warms_btw_lats  =  50000
trajecs       =  1000
traj_between_meas = 10
initial_lat   = "l164b560.20000"
output_infile = "in.hotloop"

beta_min = float(input("beta_min = "))
beta_max = float(input("beta_max = "))
beta_delta = float(input("beta_delta = "))
h_min = float(input("h_min = "))
h_max = float(input("h_max = "))
h_delta = float(input("h_delta = "))

prompt = "initial_warms = [{}]: ".format(initial_warms)
t = input(prompt)
initial_warms = int(t) if t != '' else initial_warms 

prompt = "warms_btw_lats = [{}]: ".format(warms_btw_lats)
t = input(prompt)
warms_btw_lats = int(t) if t != '' else warms_btw_lats 

prompt = "trajecs = [{}]: ".format(trajecs)
t = input(prompt)
trajecs = int(t) if t != '' else trajecs

prompt = "traj_between_meas = [{}]: ".format(traj_between_meas)
t = input(prompt)
traj_between_meas = int(t) if t != '' else traj_between_meas

prompt = "initial lattice = [{}]: ".format(initial_lat)
t = input(prompt)
initial_lat = t if t != '' else initial_lat 

prompt = "output in.file = [{}]: ".format(output_infile)
t = input(prompt)
output_infile = t if t != '' else output_infile


eprint("beta: {} to {} by {}".format(beta_min, beta_max, beta_delta))
eprint("h   : {} to {} by {}".format(h_min, h_max, h_delta))
eprint("initial warms : {}".format(initial_warms))
eprint("warms_btw_lats: {}".format(warms_btw_lats))
eprint("initial lattice: {}".format(initial_lat))
eprint("trajecs: {}".format(trajecs))
eprint("traj_between_meas: {}".format(traj_between_meas))
eprint("output written to: {}".format(output_infile))

ofile = open(output_infile, "w")

import random
header = """prompt 0
nx 16
ny 16
nz 16
nt 4
iseed {}

""".format(random.randint(0,10000))

print(header, file=ofile)

first = """warms {}
trajecs {}
traj_between_meas {}
beta {}
h {}
steps_per_trajectory 1
qhb_steps 1
reload_serial {}
no_gauge_fix
forget
""".format(initial_warms, trajecs, traj_between_meas, 
           beta_min, h_min, initial_lat)

#print(first, file=ofile)

import numpy as np
nb = round((beta_max - beta_min)/beta_delta + 1)
nh = round((h_max - h_min)/h_delta + 1)
barray = np.linspace(beta_min, beta_max, nb)
barray = barray[::-1]
for h in np.linspace(h_min, h_max, nh):
    barray = barray[::-1]
    for b in barray:
        for let in ['A','B','C','D','E','F','G','H','I','J']:
            savelat = "l164b{:3d}h{:02d}{}".format(int(round(100*b)),
                                                   int(round(100*h)),let)
            stanza = """warms {}
trajecs {}
traj_between_meas {}
beta {:3.2f}
h {:0.2f}
steps_per_trajectory 1
qhb_steps 1
hot
no_gauge_fix
save_serial {}
""".format(warms_btw_lats, trajecs, traj_between_meas, 
           b, h, savelat)
            print(stanza, file=ofile)
