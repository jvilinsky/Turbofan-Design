import numpy as np
from Parameters import *

class Turbine:
    def __init__(self):
        self.Ut = 350
        self.Ca = 200
        self.N = 180.61 
        self.To2_To1 = 435.61

    def turbine(self):
        Wc = cpa*self.To2_To1
        dTos = Wc/cpg
        psi = (2*cpg*dTos)/self.Ut**2
        