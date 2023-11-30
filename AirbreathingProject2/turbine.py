import numpy as np
from Parameters import *
import math 

class Turbine:
    def __init__(self):
        self.Ut = 350
        #Chose a Ca of 300
        #Assume C1 is axial
        self.Ca = 300

    def turb(self, T02_T01, N, T0_comb, T0_fan,P0_comb):


        
        f_ideal = ((cpg*(1/1000) * To4) - (cpa*(1/1000) * T0_comb)) / ((43100) - (cpg*(1/1000) * To4))
        mc = ma / (BPR+1)
        mf = mc * f_ideal 
        mt = mf + mc

        dT_comp = T0_comb - T0_fan
        w_tot = mc*(cpa*(1/1000) * dT_comp) + ma*(cpa*(1/1000) * T02_T01)
        delta_T_turb = w_tot / (mt*(cpg*(1/1000)))
        Tos_exact = (psi_max*self.Ut**2)/(2*cpg)
        stagest = int(math.ceil(delta_T_turb/Tos_exact))
        Tos = dT_comp/stagest
        psi = (2*cpg*Tos)/self.Ut**2
        
        ##first turbine stage
        a3 = 0.0174533 #assume a3 in radians
        phi = self.Ca/self.Ut

        C1 = self.Ca/np.cos(a3)
        T1 = To4 - C1**2/(2*cpg)
        P1 = P0_comb*(T1/To4)**(1/m1_mt)
        rho1 = P1*100000/(R*T1)
        A1 = mt/(rho1*self.Ca)

        b3 = np.arctan(np.tan(a3)+(1/phi))
        deg_reac = (2*phi*np.tan(b3) - psi/2)/2

        ##between rotor and stator
        b2 = (psi/2 - 2* deg_reac)/(2*phi)
        a2 = np.arctan(np.tan(b2) + (1/phi))
        C2 = self.Ca/np.cos(a2) 
        Cw2 = self.Ca*np.tan(a2)

        T02_T2 = C2**2/(2*cpg)
        Ti = To4 - T02_T2
        Ti_Treal = n_loss * (C2**2/(2*cpg))
        Treal = Ti - Ti_Treal

        P01_P2 = (To4/Treal)**(1/m1_mt)
        P2 = P0_comb/P01_P2
        rho2 = P2*100000/(R*Treal)
        A2 = mt/(rho2*self.Ca)
        A2N = mt/(rho2*C2)

        T03 = To4 - Tos
        P01_P03 = (T03/To4)**(1/m1_mt)
        C3 = C1
        T3 = T03 - C3**2/(2*cpg)  
        P3 = (P01_P03)*P0_comb*(T3/T03)**(1/m1_mt)
        rho3 = P3*100000/(R*T3)
        A3 = mt/(C3*rho3)
        Cw3 = self.Ca*np.tan(a3)

        rm_stg1 = self.Ut / (2*np.pi*N)
        h_stg1 = (A3*N) / self.Ut
        rt_stg1 = rm_stg1 + h_stg1/2
        rr_stg1 = rm_stg1 - h_stg1/2 
        Ut1 = 2*np.pi* rt_stg1
        Ur1 = 2*np.pi* rr_stg1
        phi_t = self.Ca/Ut1
        phi_r = self.Ca/Ur1
        
        Cw2_h = Cw2*(rr_stg1/rm_stg1)
        Cw2_t = Cw2*(rt_stg1/rm_stg1)

        Cw3_h = Cw3*(rr_stg1/rm_stg1)
        Cw3_t = Cw3*(rt_stg1/rm_stg1)

        C2_h = np.sqrt(Cw2_h**2 + self.Ca**2)
        C2_t = np.sqrt(Cw2_t**2 + self.Ca**2)

        C3_h = np.sqrt(Cw3_h**2 + self.Ca**2)
        C3_t = np.sqrt(Cw3_t**2 + self.Ca**2)

        a2_h = np.arctan((rm_stg1/rr_stg1)*np.tan(a2))
        a2_t = np.arctan((rm_stg1/rt_stg1)*np.tan(a2))

        a3_h = np.arctan((rm_stg1/rr_stg1)*np.tan(a3))
        a3_t = np.arctan((rm_stg1/rt_stg1)*np.tan(a3))

        b2_h = np.arctan((rm_stg1/rr_stg1)*np.tan(a2) + (rr_stg1/rm_stg1)*Ur1/self.Ca)      
        b2_t = np.arctan((rm_stg1/rt_stg1)*np.tan(a2) + (rt_stg1/rm_stg1)*Ur1/self.Ca)      

        b3_h = np.arctan((rm_stg1/rr_stg1)*np.tan(a3) + (rr_stg1/rm_stg1)*Ur1/self.Ca)      
        b3_t = np.arctan((rm_stg1/rt_stg1)*np.tan(a3) + (rt_stg1/rm_stg1)*Ur1/self.Ca)  

        Al2 = np.zeros(stagest)
        Al2_h = np.zeros(stagest)
        Al2_t = np.zeros(stagest)
        Al3 = np.array([1,1,1,1])
        Al3_h = np.zeros(stagest)
        Al3_t = np.zeros(stagest)

        B2 = np.zeros(stagest)
        B2_h = np.zeros(stagest)
        B2_t = np.zeros(stagest)
        B3 = np.zeros(stagest)
        B3_h = np.zeros(stagest)
        B3_t = np.zeros(stagest)

        PR = np.zeros(stagest)
        Ar1 = np.zeros(stagest)
        Ar2 = np.zeros(stagest)
        Ar3 = np.zeros(stagest)
        h = np.zeros(stagest)

        CW2 = np.zeros(stagest)
        CW2_h = np.zeros(stagest)
        CW2_t = np.zeros(stagest)
        CW3 = np.zeros(stagest)
        CW3_h = np.zeros(stagest)
        CW3_t = np.zeros(stagest)
        c2 = np.zeros(stagest)
        c2_h = np.zeros(stagest)
        c2_t = np.zeros(stagest)
        c3_h = np.zeros(stagest)
        c3_t = np.zeros(stagest)
        u = np.zeros(stagest)
        ut = np.zeros(stagest)
        uh = np.zeros(stagest)
        
        CW2[0] = Cw2
        CW2_h[0] = Cw2_h
        CW2_t[0] = Cw2_t
        CW3[0] = Cw3
        CW3_h[0] = Cw3_h
        CW3_t[0] = Cw3_t
        c2[0] = C2
        c2_h[0] = C2_h
        c2_t[0] = C2_t
        c3_h[0] = C3_h
        c3_t[0] = C2_t
        u[0] = self.Ut
        ut[0] = Ut1
        uh[0] = Ur1

        Ar1[0] = A1
        Ar2[0] = A2
        Ar3[0] = A3
        

        Al2[0] = a2
        Al2_h[0] = a2_h
        Al2_t[0] = a2_t
        Al3[0] = a3
        Al3_h[0] = a3_h
        Al3_t[0] = a3_t
        h[0] = h_stg1
        PR[0] = P01_P03
        
        B2[0] = b2
        B2_h[0] = b2_h
        B2_t[0] = b2_t
        B3[0] = b3
        B3_h[0] = b3_h
        B3_t[0] = b3_t

        for i in range(stagest-1):
            a3 = np.deg2rad(Al3[i])
            phi = self.Ca/self.Ut

            C1 = self.Ca/np.cos(a3)
            T1 = To4 - C1**2/(2*cpg)
            P1 = P0_comb*(T1/To4)**(1/m1_mt)
            rho1 = P1*100000/(R*T1)
            A1 = mt/(rho1*self.Ca)

            b3 = np.arctan(np.tan(a3)+(1/phi))
            deg_reac = (2*phi*np.tan(b3) - psi/2)/2

            ##between rotor and stator
            b2 = (psi/2 - 2* deg_reac)/(2*phi)
            a2 = np.arctan(np.tan(b2) + (1/phi))
            C2 = self.Ca/np.cos(a2) 
            Cw2 = self.Ca*np.tan(a2)

            T02_T2 = C2**2/(2*cpg)
            Ti = To4 - T02_T2
            Ti_Treal = n_loss * (C2**2/(2*cpg))
            Treal = Ti - Ti_Treal

            P01_P2 = (To4/Treal)**(1/m1_mt)
            P2 = P0_comb/P01_P2
            rho2 = P2*100000/(R*Treal)
            A2 = mt/(rho2*self.Ca)
            A2N = mt/(rho2*C2)

            T03 = To4 - Tos
            P01_P03 = (T03/To4)**(1/m1_mt)
            C3 = C1
            T3 = T03 - C3**2/(2*cpg)  
            P3 = (P01_P03)*P0_comb*(T3/T03)**(1/m1_mt)
            rho3 = P3*100000/(R*T3)
            A3 = mt/(C3*rho3)
            Cw3 = self.Ca*np.tan(a3)

            rm_stg1 = self.Ut / (2*np.pi*N)
            h_stg1 = (A3*N) / self.Ut
            rt_stg1 = rm_stg1 + h_stg1/2
            rr_stg1 = rm_stg1 - h_stg1/2 
            Ut1 = 2*np.pi* rt_stg1
            Ur1 = 2*np.pi* rr_stg1
            phi_t = self.Ca/Ut1
            phi_r = self.Ca/Ur1
            
            Cw2_h = Cw2*(rr_stg1/rm_stg1)
            Cw2_t = Cw2*(rt_stg1/rm_stg1)

            Cw3_h = Cw3*(rr_stg1/rm_stg1)
            Cw3_t = Cw3*(rt_stg1/rm_stg1)

            C2_h = np.sqrt(Cw2_h**2 + self.Ca**2)
            C2_t = np.sqrt(Cw2_t**2 + self.Ca**2)

            C3_h = np.sqrt(Cw3_h**2 + self.Ca**2)
            C3_t = np.sqrt(Cw3_t**2 + self.Ca**2)

            a2_h = np.arctan((rm_stg1/rr_stg1)*np.tan(a2))
            a2_t = np.arctan((rm_stg1/rt_stg1)*np.tan(a2))

            a3_h = np.arctan((rm_stg1/rr_stg1)*np.tan(a3))
            a3_t = np.arctan((rm_stg1/rt_stg1)*np.tan(a3))

            b2_h = np.arctan((rm_stg1/rr_stg1)*np.tan(a2) + (rr_stg1/rm_stg1)*Ur1/self.Ca)      
            b2_t = np.arctan((rm_stg1/rt_stg1)*np.tan(a2) + (rt_stg1/rm_stg1)*Ur1/self.Ca)      

            b2_h = np.arctan((rm_stg1/rr_stg1)*np.tan(a3) + (rr_stg1/rm_stg1)*Ur1/self.Ca)      


            CW2[i+1] = Cw2
            CW2_h[i+1] = Cw2_h
            CW2_t[i+1] = Cw2_t
            CW3[i+1] = Cw3
            CW3_h[i+1] = Cw3_h
            CW3_t[i+1] = Cw3_t
            c2[i+1] = C2
            c2_h[i+1] = C2_h
            c2_t[i+1] = C2_t
            c3_h[i+1] = C3_h
            c3_t[i+1] = C2_t
            u[i+1] = self.Ut
            ut[i+1] = Ut1
            uh[i+1] = Ur1
            
            Ar1[0] = A1
            Ar2[0] = A2
            Ar3[0] = A3

            Al2[i+1] = a2
            Al2_h[i+1] = a2_h
            Al2_t[i+1] = a2_t
            Al3_h[i+1] = a3_h
            Al3_t[i+1] = a3_t
            h[i+1] = h_stg1
            PR[i+1] = P01_P03
            
            B2[i+1] = b2
            B2_h[i+1] = b2_h
            B2_t[i+1] = b2_t
            B3[i+1] = b3
            B3_h[i+1] = b3_h
            B3_t[i+1] = b3_t

    data1 = np.array([['N (rev/s)', 'Number of Compressor Stages', 'Stage 1 Annulus Inlet Area (m^2)', 'Stage 1 Annulus Rotor Area (m^2)', 'Stage 1 Annulus Stator Area (m^2)','Stage 1 Annulus Height (m)', 'Stage 2 Annulus Area (m^2)',\
                        'Stage 2 Annulus Inlet Area (m^2)', 'Stage 2 Annulus Rotor Area (m^2)', 'Stage 2 Annulus Stator Area (m^2)','Stage 2 Annulus Height (m)',\
                              'Stage 3 Annulus Inlet Area (m^2)', 'Stage 3 Annulus Rotor Area (m^2)', 'Stage 3 Annulus Stator Area (m^2)','Stage 3 Annulus Height (m)',\
                                'Stage 4 Annulus Inlet Area (m^2)', 'Stage 4 Annulus Rotor Area (m^2)', 'Stage 4 Annulus Stator Area (m^2)','Stage 4 Annulus Height (m)',\
                                      'Stage 1 Rotor Mean Line Rotor Blade Angles (rad)','Stage 1 Rotor Blade Tip Angles (rad)','Stage 1 Rotor Hub Tip Angles (rad)','Stage 1 Stator Mean Line Blade Angle (rad)','Stage 2 Beta 2 (deg)', \
                                'Stage 3 Pressure Ratio','Stage 3 Alpha 1 (deg)','Stage 3 Alpha 2 (deg)','Stage 3 Beta 1 (deg)','Stage 3 Beta 2 (deg)', \
                                    'Stage 4 Pressure Ratio','Stage 4 Alpha 1 (deg)','Stage 4 Alpha 2 (deg)','Stage 4 Beta 1 (deg)','Stage 4 Beta 2 (deg)',\
                                        'Stage 5 Pressure Ratio','Stage 5 Alpha 1 (deg)','Stage 5 Alpha 2 (deg)','Stage 5 Beta 1 (deg)','Stage 5 Beta 2 (deg)',\
                                            'Stage 6 Pressure Ratio','Stage 6 Alpha 1 (deg)','Stage 6 Alpha 2 (deg)','Stage 6 Beta 1 (deg)','Stage 6 Beta 2 (deg)',\
                                                'Stage 7 Pressure Ratio','Stage 7 Alpha 1 (deg)','Stage 7 Alpha 2 (deg)','Stage 7 Beta 1 (deg)','Stage 7 Beta 2 (deg)',\
                                                        'Stage 8 Pressure Ratio','Stage 8 Alpha 1 (deg)','Stage 8 Alpha 2 (deg)','Stage 8 Beta 1 (deg)','Stage 8 Beta 2 (deg)',\
                                                            'Stage 9 Pressure Ratio','Stage 9 Alpha 1 (deg)','Stage 9 Alpha 2 (deg)','Stage 9 Beta 1 (deg)','Stage 9 Beta 2 (deg)'],
                         [N, stagest, Ar, stagesc, AA[0], H1[0],AA[1], H1[1],AA[2], H1[2],AA[3], H1[3],AA[4], H1[4],AA[5], H1[5],AA[6], H1[6],AA[7], H1[7],AA[8], H1[8],\
                          PR[0], A1[0], A2[0], B1[0], B2[0], PR[1], A1[1], A2[1], B1[1], B2[1], PR[2], A1[2], A2[2], B1[2], B2[2],PR[3], A1[3], A2[3], B1[3], B2[3],\
                             PR[4], A1[4], A2[4], B1[4], B2[4],PR[5], A1[5], A2[5], B1[5], B2[5],PR[6], A1[6], A2[6], B1[6], B2[6],PR[7], A1[7], A2[7], B1[7], B2[7],PR[8], A1[8], A2[8], B1[8], B2[8]]])


        
        print('Calcs Complete')

        
        
P0_comb = 21.49
T0_fan = 319.723
T0_comb = 755.33
T02_T01 = 31.722
N = 180.61
turbine = Turbine()
turbine.turb(T02_T01,N,T0_comb,T0_fan,P0_comb)
