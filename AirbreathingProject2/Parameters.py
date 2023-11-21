import numpy as np

BPR = 6
CPR = 15
FPR = 1.4
ma = 280 #kg/s

##For Compressor##
polyE_c = 0.9
polyE_f = 0.92
v2_v1_max = 0.65

##For Turbine##
To4 = 2275 #K
polyE_t = 0.92
phi_min = 0.78 
psi_max = 3.3

Ta = 288 #K
g_a = 1.4
g_g = 1.33
cpa = 1005 #j/kg*K
cpg = 1148 #j/kg*K
Pa = 1.01 #bar
g_g1a = g_a/(g_a-1)
R = 287 #j/kg*K
n1_nf = (1/polyE_f)*((g_a-1)/(g_a))
n1_nc = (1/polyE_c)*((g_a-1)/(g_a))

