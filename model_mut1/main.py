#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sun Aug 21 11:06:30 2022

@author: nicolaskylilis
"""

# packages
import simulation
import cell_size
import pandas as pd
import numpy as np

#%% run simulation

sim, df_new = simulation.simulate(ns=5)


#%% analysis

t = sim.t

a       = sim.y[0,:]
mr      = sim.y[1,:]
mc      = sim.y[2,:]
mq      = sim.y[3,:]
icr     = sim.y[4,:]
icc     = sim.y[5,:]
icq     = sim.y[6,:]
rmr     = sim.y[7,:]
rmc     = sim.y[8,:]
rmq     = sim.y[9,:]
r       = sim.y[10,:]
em      = sim.y[11,:]
q       = sim.y[12,:]
ribo    = sim.y[13,:]
m_ribo  = sim.y[14,:]
# ms      = sim.y[15,:]
# ics     = sim.y[16,:]
# rms     = sim.y[17,:]
# ps      = sim.y[18,:]

# cell parameters
fpath_params = "cell_paramameters.csv"
parms = pd.read_csv(fpath_params, index_col=("Parameter"))["Value"]
gmax=       parms["gmax"]
Kgamma=     parms["Kgamma"]
lenR =      parms["lenR"]
lenO =      parms["lenO"]
lenC =      parms["lenC"]
lenRibo=    parms["lenRibo"]

## cell size 
state_var = [a[-1]]
cell_mass = cell_size.model(parms,state_var)

## biosynthesis
# growth rate @steady_state
ttrate= (rmr + rmc + rmq)*(gmax*a/(Kgamma + a))
lam= (ttrate/cell_mass)

# protein mass in aa
rp_mass = ((ribo + icr + icc +icq + rmr + rmc +rmq)*lenRibo) + (r * lenR)
em_mass = em * lenC
q_mass  = q * lenO

growth_rate = lam[-1]*60
doubl_time =  (np.log(2)/lam[-1])

print("Biosynthesis rate (aa/min): " + "{:.2e}".format(ttrate[-1]/60))
print("Cell size (aa): " + "{:.2e}".format(cell_mass))
print("Growth rate (h-1): " + str(round(growth_rate,2)))
print("Doubling time (minutes): " + str(round(doubl_time,2)))

print("Intracellular energy (molecules):" + str(round(a[-1],2)) + "\n")

# save to file
with open("analysis_results.txt", "w") as f:
    f.write("---------------------\n")
    f.write("Phenotype\n")
    f.write("---------------------\n")
    f.write("Biosynthesis rate (aa/min): " + "{:.2e}".format(ttrate[-1]/60) + "\n")
    f.write("Cell size (aa): " + "{:.2e}".format(cell_mass) + "\n")
    f.write("Growth rate (h-1): " + str(round(growth_rate,2)) + "\n")
    f.write("Doubling time (minutes): " + str(round(doubl_time,2)) + "\n")
    
    f.write("\n" +"---------------------\n")
    f.write("System state\n")
    f.write("---------------------\n")
    f.write("Intracellular energy (molecules): " + str(round(a[-1],2)) + "\n")
    















