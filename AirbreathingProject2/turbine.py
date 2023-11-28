import numpy as np
from Parameters import *

class Turbine:
    def __init__(self):
        self.Ut = 350
        self.Ca = 200

    def turb(self, T02_T01, N, T0_comb, T0_fan):
        phi = self.Ca/self.Ut  #Flow coefficient
        f_ideal = ((cpg*(1/1000) * To4) - (cpa*(1/1000) * T0_comb)) / ((43100) - (cpg*(1/1000) * To4))
        mc = ma / (BPR+1)
        mf = mc * f_ideal 
        mt = mf + mc

        dT_comp = T0_comb - T0_fan
        w_tot = (mc * cpa*(1/1000) * dT_comp) + (ma * cpa*(1/1000) * T02_T01)
        delta_T_turb = w_tot / (mt * cpg*(1/1000))
        

T0_fan = 319.723
T0_comb = 755.33
T02_T01 = 435.61 
N = 180.61
turbine = Turbine()
turbine.turb(T02_T01,N,T0_comb,T0_fan)
        