import numpy as np
from Parameters import *
import math 
import plotly.graph_objects as go

class compressor:
    def __init__(self):
        self.phi_min = phi_min
        self.psi_max = psi_max
        self.Ut = 350 #m/s
        self.Ca = 200 #m/s
        self.rh_rt1 = 0.3
        self.rh_rt2 = 0.7
        self.Tos_stage1 = 50


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
        A_2 = ma/(rho2*self.Ca)
        rh1 = self.rh_rt1*rt1
        rm1 = (rt1+rh1)/2
        h_1 = A_2/(2*np.pi*rm1)
        #print('h1=',h1)
    
        
        ##Compressor##
        mc = ma / (BPR+1)
        rt2 = np.sqrt(mc/(np.pi*rho2*self.Ca*(1-self.rh_rt2**2)))
        N2 = self.Ut/(2*np.pi*rt2)
        U2t = 2*np.pi*rt2*N2
        V2t = np.sqrt(U2t**2+self.Ca**2)
        a2 = np.sqrt(g_a*R*T2)
        M2t = V2t/a2
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
        
        #(dont need these)
        '''
        ##Fan Stage Calculatioin##
        dTo21 = To2 - To1
        U1 = 2*np.pi*rm1*N1
        B1 = np.arctan(U1/self.Ca)
        V1 = self.Ca/np.cos(B1)
        V2 = V1 * v2_v1_max
        B2 = np.arccos(self.Ca/V2)
        Tos_f = (U1*self.Ca*(np.tan(B1)-np.tan(B2)))/cpa
        stagesf = dTo21/Tos_f
        '''

        ##Compressor Stage Calulation##
        dTo32 = To3 - To2
        U2 = 2*np.pi*rm2*N2
        B1 = np.arctan(U2/self.Ca)
        V1 = self.Ca/np.cos(B1)
        V2 = V1 * v2_v1_max
        B2 = np.arccos(self.Ca/V2)
        #Tos_c = (U2*self.Ca*(np.tan(B1)-np.tan(B2)))/cpa
        stagesc = int(math.ceil(dTo32/self.Tos_stage1))

        ##Stage by stage design##
        ##Fan Stage## (dont need)
        '''
        U1 = 2*np.pi*rm1*N1
        dCw = (cpa*Tos_f)/(self.lamb1*U1)
        b1 = np.arctan(U1/self.Ca)
        Cw1 = 0 
        Cw2 = dCw
        b2 = np.arctan((U1-Cw2)/self.Ca)
        a2 = np.arctan(Cw2/self.Ca)
        deHaller_rotor1 = np.cos(B1)/np.cos(B2)
        Po3_Po1 = (1+((polyE_f*Tos_f)/To1))**g_g1a
        Po3_1 = Po1*Po3_Po1
        To3_1 = To1 + Tos_f
        DOR1 = 1 - (Cw2+Cw1)/(2*U1)
        '''


        ##Stage 2 (Compressor)##
        self.Tos_stage2 = dTo32/stagesc
        dCw = (cpa*self.Tos_stage2)/(lamb2*U2)
        b1 = np.arctan(U2/self.Ca)
        Cw1 = 0
        Cw2 = dCw + Cw1
        b2 = np.arctan((U2-Cw2)/self.Ca)
        a1 = np.arctan(Cw1/self.Ca)
        a2 = np.arctan(Cw2/self.Ca)
        deHaller_rotor2 = np.cos(b1)/np.cos(b2)
        Po3_Po2 = (1+((polyE_c*self.Tos_stage2)/To2))**g_g1a
        Po3_2 = Po2*Po3_Po2
        To3_2 = To2 + self.Tos_stage2
        rt_2 = rm2 + (h2/2) # tip
        rr_2 = rm2 - (h2/2) #root
        Ut_2 = 2*np.pi*rt_2*N2
        Ur_2 = 2*np.pi*rr_2*N2
        Cw_1_t_2 = Cw1 * (rm2 / rt_2)
        Cw_1_r_2 = Cw1 * (rm2 / rr_2)
        Cw_2_t_2 = Cw2 * (rm2 / rt_2)
        Cw_2_r_2 = Cw2 * (rm2 / rr_2)
        alpha_1_t_2 = np.arctan(Cw1 / self.Ca)
        alpha_1_r_2 = np.arctan(Cw1 / self.Ca)
        alpha_2_t_2 = np.arctan(Cw2 / self.Ca)
        alpha_2_r_2 = np.arctan(Cw2 / self.Ca)
        Beta_1_t_2 = (Ut_2 - Cw1) / self.Ca
        Beta_1_r_2 = (Ur_2 - Cw1) / self.Ca
        Beta_2_t_2 = (Ut_2 - Cw2) / self.Ca
        Beta_2_r_2 = (Ur_2 - Cw2) / self.Ca
        DOR_t_2 = 1 - (Cw2 / (2 * Ut_2))
        DOR_r_2 = 1 - (Cw2 / (2 * Ur_2))
        DOR2_mean = 1 - (Cw2+Cw1)/(2*U2)

        ##Stages 3-9 (Compressor)##
        DOR = np.array([0.65,0.6,0.5,0.5,0.5,0.5,0.5,0.5])
        lamb = np.array([0.95,0.94,0.92,0.9,0.89,0.88,0.87,0.86])
        B1 = np.zeros(stagesc)
        B2 = np.zeros(stagesc)
        A1 = np.zeros(stagesc)
        A2 = np.zeros(stagesc)
        De_Haller_stator = np.zeros(8)
        De_Haller_rotor = np.zeros(8)
        TO3 = np.zeros(9)
        PO3 = np.zeros(9)
        AA = np.zeros(stagesc)
        PR = np.zeros(stagesc)
        B1[0] = b1
        B2[0] = b2
        A1[0] = a1
        A2[0] = a2
        TO3[0] = To3_2
        PO3[0] = Po3_2
        AA[0] = A3
        PR[0] = Po3_Po2

        #For Root and tip values
        H1 = np.zeros(stagesc)
        RT = np.zeros(stagesc)
        RR = np.zeros(stagesc)
        UT = np.zeros(stagesc)
        UR = np.zeros(stagesc)
        CN = np.zeros(stagesc-1)
        CW1_r = np.zeros(stagesc)
        CW1_t = np.zeros(stagesc)
        CW2_r = np.zeros(stagesc)
        CW2_t = np.zeros(stagesc)
        A1_r = np.zeros(stagesc)
        A1_t = np.zeros(stagesc)
        A2_r = np.zeros(stagesc)
        A2_t = np.zeros(stagesc)
        B1_r = np.zeros(stagesc)
        B1_t = np.zeros(stagesc)
        B2_r = np.zeros(stagesc)
        B2_t = np.zeros(stagesc)
        DOR_r = np.zeros(stagesc)
        DOR_t = np.zeros(stagesc)

        H1[0] = h2
        RT[0] = rt_2
        RR[0] = rr_2
        UT[0] = Ut_2
        UR[0] = Ur_2
        CW1_r[0] = Cw_1_r_2
        CW1_t[0] = Cw_1_t_2
        CW2_r[0] = Cw_2_r_2
        CW2_t[0] = Cw_2_t_2
        A1_r[0] = alpha_1_r_2
        A1_t[0] = alpha_1_t_2
        A2_r[0] = alpha_2_r_2
        A2_t[0] = alpha_2_t_2
        B1_r[0] = Beta_1_r_2
        B1_t[0] = Beta_1_t_2
        B2_r[0] = Beta_2_r_2
        B2_t[0] = Beta_2_t_2
        DOR_r[0] = DOR_r_2
        DOR_t[0] = DOR_t_2


        for i in range(stagesc-1):
            #Mean radius angles#
            b2 = np.arctan((-1/2)*(((self.Tos_stage2*cpa)/(lamb[i]*U2*self.Ca))-((2*DOR[i]*U2)/self.Ca)))
            b1 = np.arctan((DOR[i]*2*U2)/self.Ca-np.tan(b2))
            
            a1 = np.arctan(U2/self.Ca-np.tan(b1))
            a2 = np.arctan(U2/self.Ca-np.tan(b2))

            cw1 = self.Ca*np.tan(a1)
            cw2 = self.Ca*np.tan(a2)

            #de_Haller Values
            de_Haller_stator = np.cos(a2)/np.cos(a1)
            de_Haller_rotor = np.cos(b1)/np.cos(b2)

            Po3_Po1 = (1+(polyE_c*self.Tos_stage2)/TO3[i])**3.5
            po3 = PO3[i]*Po3_Po1
            to3 = TO3[i]+self.Tos_stage2

            #For root and tip calculations
            cn = self.Ca/np.cos(a1)
            t1 = to3 - cn**2 / (2*cpa*10**3)
            p1 = po3 * (t1/to3)**g_g1a
            rho1 = (p1*10**5)/(R*t1)

            area = mc / (rho1 * self.Ca)
            h1 = area/(2*np.pi*rm2)

            rt_stator = rm2 + (h1/2)
            rr_stator = rm2 - (h1/2)
            rt = (rt_stator+RT[i])/2
            rr = (rr_stator+RR[i])/2
            ut = 2 * np.pi*rt*N2
            ur = 2 * np.pi*rr*N2
            cw1_t = cw1 * (rm2/rt)
            cw1_r = cw1 * (rm2 / rr)
            cw2_t = cw2 * (rm2 / rt)
            cw2_r = cw2 * (rm2/ rr)
            a1_t = np.arctan(cw1_t/self.Ca)
            a1_r = np.arctan(cw1_r/self.Ca)
            a2_t = np.arctan(cw2_t/self.Ca)
            a2_r = np.arctan(cw2_r/self.Ca)
            b1_t = (ut - cw1_t)/self.Ca
            b1_r = (ur - cw1_r)/self.Ca
            b2_t = (ut - cw2_t)/self.Ca
            b2_r = (ur - cw2_r)/self.Ca
            dor_t = 1 - (cw2_t/(2*ut))
            dor_r = 1 - (cw2_r/(2*ur))

            #Filling in arrays
            B1[i+1] = b1 ; B2[i+1] = b2 ; A1[i+1] = a1 ; A2[i+1] = a2
            De_Haller_stator[i] = de_Haller_stator
            De_Haller_rotor[i] = de_Haller_rotor
            TO3[i+1] = to3 ; PO3[i+1] = po3; PR[i+1] = Po3_Po1
            H1[i+1] = h1
            AA[i+1] = area
            RT[i+1] = rt
            RR[i+1] = rr
            UT[i+1] = ut
            UR[i+1] = ur
            CN[i] = cn
            CW1_r[i+1] = cw1_r
            CW1_t[i+1] = cw1_t
            CW2_r[i+1] = cw2_r
            CW2_t[i+1] = cw2_t
            A1_r[i+1] = a1_r
            A1_t[i+1] = a1_t
            A2_r[i+1] = a2_r
            A2_t[i+1] = a2_t
            B1_r[i+1] = b1_r
            B1_t[i+1] = b1_t
            B2_r[i+1] = b2_r
            B2_t[i+1] = b2_t
            DOR_r[i+1] = dor_r
            DOR_t[i+1] = dor_t
        
        data1 = np.array([['N (rev/s)','Annulus Area (m^2)','Annulus Height (m)', 'Number of Compressor Stages', \
                           'Stage 1 Annulus Area (m^2)', 'Stage 1 Annulus Height (m)', 'Stage 2 Annulus Area (m^2)', 'Stage 2 Annulus Height (m)',\
                           'Stage 3 Annulus Area (m^2)', 'Stage 3 Annulus Height (m)','Stage 4 Annulus Area (m^2)', 'Stage 4 Annulus Height (m)',\
                           'Stage 5 Annulus Area (m^2)', 'Stage 5 Annulus Height (m)','Stage 6 Annulus Area (m^2)', 'Stage 6 Annulus Height (m)',\
                           'Stage 7 Annulus Area (m^2)', 'Stage 7 Annulus Height (m)','Stage 8 Annulus Area (m^2)', 'Stage 8 Annulus Height (m)','Stage 9 Annulus Area (m^2)', 'Stage 9 Annulus Height (m)',\
                          'Stage 1 Pressure Ratio','Stage 1 Alpha 1 (deg)','Stage 1 Alpha 2 (deg)','Stage 1 Beta 1 (deg)','Stage 1 Beta 2 (deg)', \
                            'Stage 2 Pressure Ratio','Stage 2 Alpha 1 (deg)','Stage 2 Alpha 2 (deg)','Stage 2 Beta 1 (deg)','Stage 2 Beta 2 (deg)', \
                                'Stage 3 Pressure Ratio','Stage 3 Alpha 1 (deg)','Stage 3 Alpha 2 (deg)','Stage 3 Beta 1 (deg)','Stage 3 Beta 2 (deg)', \
                                    'Stage 4 Pressure Ratio','Stage 4 Alpha 1 (deg)','Stage 4 Alpha 2 (deg)','Stage 4 Beta 1 (deg)','Stage 4 Beta 2 (deg)',\
                                        'Stage 5 Pressure Ratio','Stage 5 Alpha 1 (deg)','Stage 5 Alpha 2 (deg)','Stage 5 Beta 1 (deg)','Stage 5 Beta 2 (deg)',\
                                            'Stage 6 Pressure Ratio','Stage 6 Alpha 1 (deg)','Stage 6 Alpha 2 (deg)','Stage 6 Beta 1 (deg)','Stage 6 Beta 2 (deg)',\
                                                'Stage 7 Pressure Ratio','Stage 7 Alpha 1 (deg)','Stage 7 Alpha 2 (deg)','Stage 7 Beta 1 (deg)','Stage 7 Beta 2 (deg)',\
                                                        'Stage 8 Pressure Ratio','Stage 8 Alpha 1 (deg)','Stage 8 Alpha 2 (deg)','Stage 8 Beta 1 (deg)','Stage 8 Beta 2 (deg)',\
                                                            'Stage 9 Pressure Ratio','Stage 9 Alpha 1 (deg)','Stage 9 Alpha 2 (deg)','Stage 9 Beta 1 (deg)','Stage 9 Beta 2 (deg)'],
                         [N2, A_2, h_1, stagesc, AA[0], H1[0],AA[1], H1[1],AA[2], H1[2],AA[3], H1[3],AA[4], H1[4],AA[5], H1[5],AA[6], H1[6],AA[7], H1[7],AA[8], H1[8],\
                          PR[0], A1[0], A2[0], B1[0], B2[0], PR[1], A1[1], A2[1], B1[1], B2[1], PR[2], A1[2], A2[2], B1[2], B2[2],PR[3], A1[3], A2[3], B1[3], B2[3],\
                             PR[4], A1[4], A2[4], B1[4], B2[4],PR[5], A1[5], A2[5], B1[5], B2[5],PR[6], A1[6], A2[6], B1[6], B2[6],PR[7], A1[7], A2[7], B1[7], B2[7],PR[8], A1[8], A2[8], B1[8], B2[8]]])

        fig = go.Figure(data=[go.Table(
            header=dict(values=['Variable', 'Value'],
                        line_color='darkslategray',
                        fill_color='lightskyblue',
                        align='left'),
            cells=dict(values=data1, 
                    line_color='darkslategray',
                    fill_color='lightcyan',
                    align='left'))])
        fig.update_layout(title='Compressor Annulus and Mean Radius Calculations')
        fig.update_layout(width=800, height=1000)
        fig.show()

        data2 = np.array([['Stage 1 U tip (m/s)','Stage 1 U root (m/s)','Stage 1 Cw1 tip (m/s)','Stage 1 Cw1 root (m/s)','Stage 1 Cw2 tip (m/s)','Stage 1 Cw2 root (m/s)',\
                           "Stage 1 Alpha 1 tip (deg)",'Stage 1 Alpha 1 root (deg)','Stage 1 Alpha 2 tip (deg)','Stage 1 Alpha 2 root','Stage 1 Beta 1 tip (deg)','Stage 1 Beta 1 root (deg)', 'Stage 1 Beta 2 tip (deg)','Stage 1 Beta 2 tip (deg)',\
                            'Stage 2 U tip (m/s)','Stage 2 U root (m/s)','Stage 2 Cw1 tip (m/s)','Stage 2 Cw1 root (m/s)','Stage 2 Cw2 tip (m/s)','Stage 2 Cw2 root (m/s)',\
                           "Stage 2 Alpha 1 tip (deg)",'Stage 2 Alpha 1 root (deg)','Stage 2 Alpha 2 tip (deg)','Stage 2 Alpha 2 root','Stage 2 Beta 1 tip (deg)','Stage 2 Beta 1 root (deg)', 'Stage 2 Beta 2 tip (deg)','Stage 2 Beta 2 tip (deg)',\
                            'Stage 3 U tip (m/s)','Stage 3 U root (m/s)','Stage 3 Cw1 tip (m/s)','Stage 3 Cw1 root (m/s)','Stage 3 Cw2 tip (m/s)','Stage 3 Cw2 root (m/s)',\
                           "Stage 3 Alpha 1 tip (deg)",'Stage 3 Alpha 1 root (deg)','Stage 3 Alpha 2 tip (deg)','Stage 3 Alpha 2 root','Stage 3 Beta 1 tip (deg)','Stage 3 Beta 1 root (deg)', 'Stage 3 Beta 2 tip (deg)','Stage 3 Beta 2 tip (deg)',\
                            'Stage 4 U tip (m/s)','Stage 4 U root (m/s)','Stage 4 Cw1 tip (m/s)','Stage 4 Cw1 root (m/s)','Stage 4 Cw2 tip (m/s)','Stage 4 Cw2 root (m/s)',\
                           "Stage 4 Alpha 1 tip (deg)",'Stage 4 Alpha 1 root (deg)','Stage 4 Alpha 2 tip (deg)','Stage 4 Alpha 2 root','Stage 4 Beta 1 tip (deg)','Stage 4 Beta 1 root (deg)', 'Stage 4 Beta 2 tip (deg)','Stage 4 Beta 2 tip (deg)',\
                            'Stage 5 U tip (m/s)','Stage 5 U root (m/s)','Stage 5 Cw1 tip (m/s)','Stage 5 Cw1 root (m/s)','Stage 5 Cw2 tip (m/s)','Stage 5 Cw2 root (m/s)',\
                           "Stage 5 Alpha 1 tip (deg)",'Stage 5 Alpha 1 root (deg)','Stage 5 Alpha 2 tip (deg)','Stage 5 Alpha 2 root','Stage 5 Beta 1 tip (deg)','Stage 5 Beta 1 root (deg)', 'Stage 5 Beta 2 tip (deg)','Stage 5 Beta 2 tip (deg)',\
                            'Stage 6 U tip (m/s)','Stage 6 U root (m/s)','Stage 6 Cw1 tip (m/s)','Stage 6 Cw1 root (m/s)','Stage 6 Cw2 tip (m/s)','Stage 6 Cw2 root (m/s)',\
                           "Stage 6 Alpha 1 tip (deg)",'Stage 6 Alpha 1 root (deg)','Stage 6 Alpha 2 tip (deg)','Stage 6 Alpha 2 root','Stage 6 Beta 1 tip (deg)','Stage 6 Beta 1 root (deg)', 'Stage 6 Beta 2 tip (deg)','Stage 6 Beta 2 tip (deg)',\
                            'Stage 7 U tip (m/s)','Stage 7 U root (m/s)','Stage 7 Cw1 tip (m/s)','Stage 7 Cw1 root (m/s)','Stage 7 Cw2 tip (m/s)','Stage 7 Cw2 root (m/s)',\
                           "Stage 7 Alpha 1 tip (deg)",'Stage 7 Alpha 1 root (deg)','Stage 7 Alpha 2 tip (deg)','Stage 7 Alpha 2 root','Stage 7 Beta 1 tip (deg)','Stage 7 Beta 1 root (deg)', 'Stage 7 Beta 2 tip (deg)','Stage 7 Beta 2 tip (deg)',\
                            'Stage 8 U tip (m/s)','Stage 8 U root (m/s)','Stage 8 Cw1 tip (m/s)','Stage 8 Cw1 root (m/s)','Stage 8 Cw2 tip (m/s)','Stage 8 Cw2 root (m/s)',\
                           "Stage 8 Alpha 1 tip (deg)",'Stage 8 Alpha 1 root (deg)','Stage 8 Alpha 2 tip (deg)','Stage 8 Alpha 2 root','Stage 8 Beta 1 tip (deg)','Stage 8 Beta 1 root (deg)', 'Stage 8 Beta 2 tip (deg)','Stage 8 Beta 2 tip (deg)',\
                            'Stage 9 U tip (m/s)','Stage 9 U root (m/s)','Stage 9 Cw1 tip (m/s)','Stage 9 Cw1 root (m/s)','Stage 9 Cw2 tip (m/s)','Stage 9 Cw2 root (m/s)',\
                           "Stage 9 Alpha 1 tip (deg)",'Stage 9 Alpha 1 root (deg)','Stage 9 Alpha 2 tip (deg)','Stage 9 Alpha 2 root','Stage 9 Beta 1 tip (deg)','Stage 9 Beta 1 root (deg)', 'Stage 9 Beta 2 tip (deg)','Stage 9 Beta 2 tip (deg)'],
                          [UT[0], UR[0], CW1_t[0], CW1_r[0], CW2_t[0], CW2_r[0], A1_t[0], A1_r[0], A2_t[0], A2_r[0], B1_t[0], B1_r[0], B2_t[0], B2_r[0],\
                           UT[1], UR[1], CW1_t[1], CW1_r[1], CW2_t[1], CW2_r[1], A1_t[1], A1_r[1], A2_t[1], A2_r[1], B1_t[1], B1_r[1], B2_t[1], B2_r[1],\
                           UT[2], UR[2], CW1_t[2], CW1_r[2], CW2_t[2], CW2_r[2], A1_t[2], A1_r[2], A2_t[2], A2_r[2], B1_t[2], B1_r[2], B2_t[2], B2_r[2],\
                           UT[3], UR[3], CW1_t[3], CW1_r[3], CW2_t[3], CW2_r[3], A1_t[3], A1_r[3], A2_t[3], A2_r[3], B1_t[3], B1_r[3], B2_t[3], B2_r[3],\
                           UT[4], UR[4], CW1_t[4], CW1_r[4], CW2_t[4], CW2_r[4], A1_t[4], A1_r[4], A2_t[4], A2_r[4], B1_t[4], B1_r[4], B2_t[4], B2_r[4],\
                           UT[5], UR[5], CW1_t[5], CW1_r[5], CW2_t[5], CW2_r[5], A1_t[5], A1_r[5], A2_t[5], A2_r[5], B1_t[5], B1_r[5], B2_t[5], B2_r[5],\
                           UT[6], UR[6], CW1_t[6], CW1_r[6], CW2_t[6], CW2_r[6], A1_t[6], A1_r[6], A2_t[6], A2_r[6], B1_t[6], B1_r[6], B2_t[6], B2_r[6],\
                           UT[7], UR[7], CW1_t[7], CW1_r[7], CW2_t[7], CW2_r[7], A1_t[7], A1_r[7], A2_t[7], A2_r[7], B1_t[7], B1_r[7], B2_t[7], B2_r[7],\
                           UT[8], UR[8], CW1_t[8], CW1_r[8], CW2_t[8], CW2_r[8], A1_t[8], A1_r[8], A2_t[8], A2_r[8], B1_t[8], B1_r[8], B2_t[8], B2_r[8]]])
        
        fig = go.Figure(data=[go.Table(
            header=dict(values=['Variable', 'Value'],
                        line_color='darkslategray',
                        fill_color='lightskyblue',
                        align='left'),
            cells=dict(values=data2, 
                    line_color='darkslategray',
                    fill_color='lightcyan',
                    align='left'))])
        fig.update_layout(title='Compressor Tip and Root Calculations')
        fig.update_layout(width=800, height=1000)
        fig.show()

Compressor = compressor()
Compressor.comp()

