import numpy as np
from Parameters import *

class Compression:

    def __init__(self):
        self.BPR = BPR

    def Totals(self):
        Ut = 350 # m/s
        Ca = 200 # m/s   assume Ca = constant
        rh_rt_f = .3  # fan
        rh_rt_c = .78  # compressor

        # before the fan
        T1 = T01 - Ca**2 / (2 * cpa * 10**3)
        p1 = P01 * (T1/T01)**y_y_1a
        rho1 = (p1 * 100) / (R * T1)

        # after the fan
        P02 = PR_fan * P01
        T02 = T01 * PR_fan**n_1_nf
        T2 = T02 - Ca**2/(2 * cpa * 10**3)
        p2 = P02 * (T2/T02)**y_y_1a
        rho2 = (p2 * 100)/(R*T2)

        # beginning of compressor
        P03 = PR_cum * P02
        T03 = T02 * PR_cum ** n_1_nc
        T3 = T03 - Ca ** 2 / (2 * cpa * 10 ** 3)
        p3 = P03 * (T3 / T03) ** y_y_1a
        rho3 = (p3 * 100) / (R * T3)
        rt3 = np.sqrt(mdot_c / (np.pi * rho2 * Ca * (1 - rh_rt_c ** 2)))

        N3 = Ut / (2 * np.pi * rt3)
        Ut3 = 2 * np.pi * rt3 * N3
        V3t = np.sqrt(Ut3 ** 2 + Ca ** 2)
        a3 = np.sqrt(ya * R * 1000 * T2)
        M3t = V3t / a3
        print(M3t)

        # dimensions of compressor stg last
        A3 = mdot_c / (rho2 * Ca)
        rh_c = rh_rt_c * rt3
        rm_c = (rh_c + rt3) / 2   # constant throughout
        h3 = A3 / (2 * np.pi * rm_c)
        rt_3 = rm_c + (h3 / 2)  # tip
        rr_3 = rm_c - (h3 / 2)  # root

        T02_T01 = T02 - T01
        T03_T02 = T03 - T02
        Um_c = 2 * np.pi * rm_c * N3
        Ur_c = 2 * np.pi * rr_3 * N3

        return T02_T01, T03_T02, Um_c, Ca, T02, P02, rh_rt_c, Ut

    def numberstages(self, T02_T01, T03_T02, Um_c, Ca):

        Bc1 = np.arctan(Um_c / Ca)
        Vc1 = Ca/np.cos(Bc1)
        Vc2 = Vc1 * de_Haller
        Bc2 = np.arccos(Ca/Vc2)
        # delta_T_c_per = Um_c * Ca * (np.tan(Bc1) - np.tan(Bc2)) / (cpa * 10**3)
        delta_T_c_per = 50    # guess
        stages_com = T03_T02 / delta_T_c_per
        print(stages_com, delta_T_c_per)
        stages_com = int(np.round(stages_com))
        delta_T_c_per = T03_T02 / stages_com

        return delta_T_c_per, stages_com

    def stage(self, delta_T_c_per, Um_c, Ca, stages_com, T02, P02, rh_rt_c, Ut):

        # stage 1
        lambda1 = .98
        Cw2 = 0
        delta_T_c_per_1 = delta_T_c_per + 5
        delta_T_c_per = (stages_com * delta_T_c_per - delta_T_c_per_1)/ (stages_com - 1)

        delta_C_w_c = (cpa * 10 ** 3 * delta_T_c_per_1) / (lambda1 * Um_c)
        Cw3 = delta_C_w_c + Cw2
        Bc1 = np.arctan(Um_c / Ca)
        Bc2 = np.arctan((Um_c - Cw3) / Ca)
        alpha1 = np.arctan((Um_c/Ca) - np.tan(Bc1))
        alpha2 = np.arctan((Um_c/Ca) - np.tan(Bc2))
        Vc3_Vc2 = np.cos(Bc1) / np.cos(Bc2)
        P03_P02 = (1 + (n_inf_c * delta_T_c_per_1 / T02)) ** y_y_1a
        P03 = P02 * P03_P02
        T03 = T02 + delta_T_c_per_1
        deg_react_f2 = 1 - (Cw2 + Cw3) / (2 * Um_c)
        print('Com stage 1 deg of react', Vc3_Vc2)


        # root and tip calcs
        T3 = T03 - Ca ** 2 / (2 * cpa * 10 ** 3)
        p3 = P03 * (T3 / T03) ** y_y_1a
        rho3 = (p3 * 100) / (R * T3)
        r2t = np.sqrt(mdot_c / (np.pi * rho3 * Ca * (1 - rh_rt_c ** 2)))
        N3 = Ut / (2 * np.pi * r2t)  # Ut is tip U from previous iteration

        # dimensions for inlet compressor stg 1
        A_1 = mdot_c / (rho3 * Ca)
        rh_c_1 = rh_rt_c * r2t
        rm_c_1 = (rh_c_1 + r2t) / 2
        h_1 = A_1 / (2 * np.pi * rm_c_1)
        rt_1 = rm_c_1 + (h_1 / 2)  # tip
        rr_1 = rm_c_1 - (h_1 / 2)  # root
        Ut_1 = 2 * np.pi * rt_1 * N3  # finding new Ut
        Um_1 = 2 * np.pi * rm_c_1 * N3
        Ur_1 = 2 * np.pi * rr_1 * N3

        Cw_1_t_1 = C2w2 * (rm_c_1 / rt_1)
        Cw_1_r_1 = Cw2 * (rm_c_1 / rr_1)
        Cw_2_t_1 = Cw3 * (rm_c_1 / rt_1)
        Cw_2_r_1 = Cw3 * (rm_c_1 / rr_1)
        alpha_1_t_1 = np.arctan(Cw2 / Ca)
        alpha_1_r_1 = np.arctan(Cw2 / Ca)
        alpha_2_t_1 = np.arctan(Cw3 / Ca)
        alpha_2_r_1 = np.arctan(Cw3 / Ca)
        Beta_1_t_1 = (Ut_1 - Cw2) / Ca
        Beta_1_r_1 = (Ur_1 - Cw2) / Ca
        Beta_2_t_1 = (Ut_1 - Cw3) / Ca
        Beta_2_r_1 = (Ur_1 - Cw3) / Ca
        deg_react_t_1 = 1 - (Cw3 / (2 * Ut_1))
        deg_react_r_1 = 1 - (Cw3 / (2 * Ur_1))
        # after stator
        '''C_s = Ca / np.cos(alpha1)
        T_s = T03 - (C_s**2/(2 * cpa * 10**3))
        p_s = P03*(T_s/T03)**y_y_1a
        rho_s = (p_s * 100) / (R * T_s)
        A_s = mdot_c / (rho_s * Ca)
        h_s = A_s / (2 * np.pi * rm_c_1)

        rt_s = rm_c_1 + (h_s / 2)  # tip
        rr_s = rm_c_1 - (h_s / 2)  # root
        # after rotor
        rt_2 = (rt_1 + rt_s)/2
        rr_2 = (rr_1 + rr_s)/2
        Ut_2 = 2 * np.pi * rt_2 * N3
        Ur_2 = 2 * np.pi * rr_2 * N3'''


        T_1_e = np.zeros(stages_com)
        p_1_e = np.zeros(stages_com)
        rho_1_e = np.zeros(stages_com)
        A_1_e = np.zeros(stages_com)
        h_1_e = np.zeros(stages_com)
        rt_1_e = np.zeros(stages_com)
        rr_1_e = np.zeros(stages_com)
        Ut_1_e = np.zeros(stages_com)
        Ur_1_e = np.zeros(stages_com)
        C_n = np.zeros(stages_com)
        Cw_1_r = np.zeros(stages_com)
        Cw_1_t = np.zeros(stages_com)
        Cw_2_r = np.zeros(stages_com)
        Cw_2_t = np.zeros(stages_com)
        alpha_1_r = np.zeros(stages_com)
        alpha_1_t = np.zeros(stages_com)
        alpha_2_r = np.zeros(stages_com)
        alpha_2_t = np.zeros(stages_com)
        Beta_1_r = np.zeros(stages_com)
        Beta_1_t = np.zeros(stages_com)
        Beta_2_r = np.zeros(stages_com)
        Beta_2_t = np.zeros(stages_com)
        deg_react_r = np.zeros(stages_com)
        deg_react_t = np.zeros(stages_com)

        deg_react_guess = np.array([.75, .7, .65, .6, .55, .53, .5, .5, .5, .5, .5])
        lambda_guess = np.array([.94, .91, .89, .87, .84, .83, .83, .83, .83, .83, .83])
        B2 = np.zeros(stages_com)
        B1 = np.zeros(stages_com)
        alpha_1 = np.zeros(stages_com)
        alpha_2 = np.zeros(stages_com)
        Cw1 = np.zeros(stages_com)
        Cw2 = np.zeros(stages_com)
        de_Haller_rotor = np.zeros(stages_com)
        de_Haller_stator = np.zeros(stages_com)
        P03_P01_2 = np.zeros(stages_com)
        P03_2_2 = np.zeros(stages_com)
        T03_2 = np.zeros(stages_com)

        for i in range(stages_com - 1):
            B1[0] = Bc1
            B2[0] = Bc2
            alpha_1[0] = alpha1
            alpha_2[0] = alpha2
            P03_P01_2[0] = P03_P02
            P03_2_2[0] = P03
            T03_2[0] = T03
            Ut_1_e[0] = Ut_1
            Ur_1_e[0] = Ur_1
            rt_1_e[0] = rt_1
            rr_1_e[0] = rr_1
            Ut_1_e[0] = Ut_1
            Ur_1_e[0] = Ur_1
            Cw_1_t[0] =Cw_1_t_1
            Cw_1_r[0] =   Cw_1_r_1
            Cw_2_t[0] = Cw_2_t_1
            Cw_2_r[0] =Cw_2_r_1
            alpha_1_t[0] = alpha_1_t_1
            alpha_1_r[0] = alpha_1_r_1
            alpha_2_t[0] = alpha_2_t_1
            alpha_2_r[0] = alpha_2_r_1
            Beta_1_t[0] = Beta_1_t_1
            Beta_1_r[0] = Beta_1_r_1
            Beta_2_t[0] = Beta_2_t_1
            Beta_2_r[0] = Beta_2_r_1
            deg_react_t[0] = deg_react_t_1
            deg_react_r[0] = deg_react_r_1


            i = i + 1
            B2[i] = np.arctan(-.5*(((delta_T_c_per * cpa * 10**3)/(lambda_guess[i] * Um_c * Ca)) - ((2 * deg_react_guess[i] * Um_c)/Ca)))
            B1[i] = np.arctan(((2*Um_c*deg_react_guess[i])/Ca) - np.tan(B2[i]))
            alpha_1[i] = np.arctan(Um_c/Ca - np.tan(B1[i]))
            alpha_2[i] = np.arctan(Um_c/Ca - np.tan(B2[i]))

            Cw1[i] = Ca * np.tan(alpha_1[i])
            Cw2[i] = Ca * np.tan(alpha_2[i])

            de_Haller_rotor[i] = np.cos(B1[i])/np.cos(B2[i])
            de_Haller_stator[i] = np.cos(alpha_2[i])/np.cos(alpha_1[i])
            P03_P01_2[i] = (1 + n_inf_c * delta_T_c_per/T03)**y_y_1a

            P03_2_2[i] = P03_2_2[i - 1] * P03_P01_2[i]
            T03_2[i] = T03_2[i - 1] + delta_T_c_per


            # root and tip calcs
            C_n[i] = Ca/np.cos(alpha_1[i])
            T_1_e[i] = T03_2[i] - C_n[i] ** 2 / (2 * cpa * 10 ** 3)
            p_1_e[i] = P03_2_2[i] * (T_1_e[i] / T03_2[i]) ** y_y_1a
            rho_1_e[i] = (p_1_e[i] * 100) / (R * T_1_e[i])

            # dimensions of compressor stg 1
            A_1_e[i] = mdot_c / (rho_1_e[i] * Ca)
            h_1_e[i] = A_1_e[i]/(2 * np.pi * rm_c_1)

            rt_stator = rm_c_1 + (h_1_e[i] / 2)
            rr_stator = rm_c_1 - (h_1_e[i] / 2)
            rt_1_e[i] = (rt_stator + rt_1_e[i - 1]) / 2
            rr_1_e[i] = (rr_stator + rr_1_e[i - 1]) / 2
            Ut_1_e[i] = 2 * np.pi * rt_1_e[i] * N3  # finding new Ut
            Ur_1_e[i] = 2 * np.pi * rr_1_e[i] * N3

            Cw_1_t[i] = Cw1[i] * (rm_c_1 / rt_1_e[i])
            Cw_1_r[i] = Cw1[i] * (rm_c_1 / rr_1_e[i])
            Cw_2_t[i] = Cw2[i] * (rm_c_1 / rt_1_e[i])
            Cw_2_r[i] = Cw2[i] * (rm_c_1 / rr_1_e[i])  # double check
            alpha_1_t[i] = np.arctan(Cw_1_t[i] / Ca)
            alpha_1_r[i] = np.arctan(Cw_1_r[i] / Ca)
            alpha_2_t[i] = np.arctan(Cw_2_t[i]/Ca)
            alpha_2_r[i] = np.arctan(Cw_2_r[i]/Ca)
            Beta_1_t[i] = (Ut_1_e[i] - Cw_1_t[i]) / Ca
            Beta_1_r[i] = (Ur_1_e[i] - Cw_1_r[i]) / Ca
            Beta_2_t[i] = (Ut_1_e[i] - Cw_2_t[i]) / Ca  # double check
            Beta_2_r[i] = (Ur_1_e[i] - Cw_2_r[i]) / Ca
            deg_react_t[i] = 1 - (Cw_2_t[i]/(2 * Ut_1_e[i]))  # double check
            deg_react_r[i] = 1 - (Cw_2_r[i]/(2 * Ur_1_e[i]))


        print(de_Haller_rotor, P03_P01_2)
