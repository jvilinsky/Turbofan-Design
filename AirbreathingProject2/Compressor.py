import numpy as np
from Parameters import *

class compressor:
    def __init__(self):
        self.phi_min = phi_min
        self.psi_max = psi_max
        self.Ut = 350 #m/s
        self.Ca = 175 #m/s
        self.rh_rt1 = 0.3
        self.rh_rt2 = 0.18
        self.lamb1 = 0.98
        self.lamb2 = 0.95


    def comp(self):
        ##Before Fan##
        To1 = Ta
        Po1 = Pa
        T1 = To1 - self.Ca**2/(2*cpa)
        P1 = Po1*(T1/To1)**g_g1a
        rho1 = (P1*10**5)/(R*T1)

        rt1 = np.sqrt(ma/(np.pi*rho1*self.Ca*(1-self.rh_rt1**2)))
        N1 = self.Ut/(2*np.pi*rt1)
        U1t = 2*np.pi*rt1*N1
        V1t = np.sqrt(U1t**2+self.Ca**2)
        a1 = np.sqrt(g_a*R*T1)
        M1t = V1t/a1

        ##Fan##
        To2 = To1*(FPR)**n1_nf
        T2 = To2 - self.Ca**2/(2*cpa)
        Po2 = FPR*Po1
        P2 = Po2*(T2/To2)**g_g1a
        rho2 = (P2*10**5)/(R*T2)
        A2 = ma/(rho2*self.Ca)
        rh1 = self.rh_rt1*rt1
        rm1 = (rt1+rh1)/2
        h1 = A2/(2*np.pi*rm1)
        #print('h1=',h1)

        mc = ma / (BPR+1)
        rt2 = np.sqrt(mc/(np.pi*rho2*self.Ca*(1-self.rh_rt2**2)))
        N2 = self.Ut/(2*np.pi*rt2)
        U2t = 2*np.pi*rt2*N2
        V2t = np.sqrt(U2t**2+self.Ca**2)
        a2 = np.sqrt(g_a*R*T2)
        M2t = V2t/a2

        ##Compressor##
        To3 = To2*(CPR)**n1_nc
        T3 = To3 - self.Ca**2/(2*cpa)
        Po3 = CPR*Po2
        P3 = Po3*(T3/To3)**g_g1a
        rho3 = (P3*10**5)/(R*T3)
        A3 = mc/(rho3*self.Ca)
        rh2 = self.rh_rt2*rt2
        rm2 = (rt2+rh2)/2
        h2 = A3/(2*np.pi*rm2)
        #print('h2=',h2)

        ##Fan Stage Calculatioin##
        dTo21 = To2 - To1
        U1 = 2*np.pi*rm1*N1
        B1 = np.arctan(U1/self.Ca)
        V1 = self.Ca/np.cos(B1)
        V2 = V1 * v2_v1_max
        B2 = np.arccos(self.Ca/V2)
        Tos_f = (U1*self.Ca*(np.tan(B1)-np.tan(B2)))/cpa
        stagesf = dTo21/Tos_f

        ##Compressor Stage Calulation##
        dTo32 = To3 - To2
        U2 = 2*np.pi*rm2*N2
        B1 = np.arctan(U2/self.Ca)
        V1 = self.Ca/np.cos(B1)
        V2 = V1 * v2_v1_max
        B2 = np.arccos(self.Ca/V2)
        Tos_c = (U2*self.Ca*(np.tan(B1)-np.tan(B2)))/cpa
        stagesc = dTo32/Tos_c

        ##Stage by stage design##
        ##Fan Stage##
        U1 = 2*np.pi*rm1*N1
        dCw = (cpa*Tos_f)/(self.lamb1*U1)
        B1 = np.arctan(U1/self.Ca)
        Cw1 = 0 
        Cw2 = dCw
        B2 = np.arctan((U1-Cw2)/self.Ca)
        a2 = np.arctan(Cw2/self.Ca)
        deHaller_rotor1 = np.cos(B1)/np.cos(B2)
        Po3_Po1 = (1+((polyE_f*Tos_f)/To1))**g_g1a
        Po3_1 = Po1*Po3_Po1
        To3_1 = To1 + Tos_f
        DOR1 = 1 - (Cw2+Cw1)/(2*U1)


        ##Stage 2 (Compressor)##
        dCw = (cpa*Tos_c)/(self.lamb2*U2)
        B1 = np.arctan(U2/self.Ca)
        Cw1 = Cw2
        Cw2 = dCw + Cw1
        B2 = np.arctan((U2-Cw2)/self.Ca)
        a2 = np.arctan(Cw2/self.Ca)
        deHaller_rotor2 = np.cos(B1)/np.cos(B2)
        Po3_Po2 = (1+((polyE_c*Tos_c)/To2))**g_g1a
        Po3_2 = Po1*Po3_Po1
        To3_2 = To2 + Tos_c
        DOR2 = 1 - (Cw2+Cw1)/(2*U1)

        ##Stages 3-9 (Compressor)##
        DOR = np.array([0.65,0.6,0.5,0.5,0.5,0.5,0.5])
        lamb = np.array([0.94,0.93,0.91,0.9,0.89,0.88,0.87])
        B1 = np.zeros(7)
        B2 = np.zeros(7)
        A1 = np.zeros(7)
        A2 = np.zeros(7)
        De_Haller_stator = np.zeros(7)
        De_Haller_rotor = np.zeros(7)
        TO3 = np.zeros(8)
        PO3 = np.zeros(8)
        TO3[0] = To3_2
        PO3[0] = Po3_2
        for i in range(7):
            b2 = np.arctan((-1/2)*(((Tos_c*cpa)/(lamb[i]*U2*self.Ca))-((2*DOR[i]*U2)/self.Ca)))
            b1 = np.arctan((DOR[i]*2*U2)/self.Ca-np.tan(b2))
            
            a1 = np.arctan(U2/self.Ca-np.tan(b1))
            a2 = np.arctan(U2/self.Ca-np.tan(b2))

            cw1 = self.Ca*np.tan(a1)
            cw2 = self.Ca*np.tan(a2)

            de_Haller_stator = np.cos(a2)/np.cos(a1)
            de_Haller_rotor = np.cos(b1)/np.cos(b2)

            Po3_Po1 = (1+(polyE_c*Tos_c)/TO3[i])**3.5
            po3 = PO3[i]*Po3_Po1
            to3 = TO3[i]+Tos_c

            B1[i] = b1 ; B2[i] = b2 ; A1[i] = a1 ; A2[i] = a2
            De_Haller_stator[i] = de_Haller_stator
            De_Haller_rotor[i] = de_Haller_rotor
            TO3[i+1] = to3
            PO3[i+1] = po3
        
        print('Calcs Completed')



Compressor = compressor()
Compressor.comp()

