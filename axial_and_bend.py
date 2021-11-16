"""
sources of knoweladge
pietryga [1]
https://chodor-projekt.net/encyclopedia/krzywe-interakcji-m-n-zelbetu/ [2]
knauf
łapko jensen
oleszek
"""

import csv
import math
import numpy as np

from t_sect_ben_reinf import TCrReinf
from char_geom import CharGeom

class GeneralAxBend(TCrReinf):
    """cross section axial force + bending moment capasity computations;
    returns stresses in conrete and steel layers according to computed itaratively
    (gradient descent) epsilon and fi (strain and rotation of cross section) 
    to achive close resulting M and N as M_sd and N_sd"""
    EC2 = 0.0020  # concrete strain of section in the begining of flat (plastic) stresses
    ECU2 = 0.0035  # max. concrete strain of section
    EYD = 0.00217  # max. steel strain (FOR BOTH B500SP AND BST500S STEEL)
    _R = 20  # curve parameter (FOR BOTH B500SP AND BST500S STEEL ?)
    E_STELL = 200 * 10 ** 3  # acc. to [2]
    EPS_INIT = 0.01 * 10 ** -3

    def __init__(self, name, b, h, hsl, 
                 beff, cl_conc, cl_steel, 
                 c, fi, fi_s, fi_opp, 
                 nl_reinf_top, nl_reinf_bottom,
                 m_sd, n_sd):
        super().__init__(name, b, h, hsl, beff, cl_conc, cl_steel, c, fi, fi_s, fi_opp, m_sd)
        self.n_sd = n_sd
        char_geom = CharGeom()
        self.e_vert = char_geom.find_center_m((b[0], h[0], b[1], h[1], b[2], h[2]))  # = sefl.h_bottom
        self.h_top = sum(self.h) - self.e_vert
        self.nl_reinf_top = nl_reinf_top
        self.nl_reinf_bottom = nl_reinf_bottom
        
    @property
    def E_cm(self):
        """returns E_cm value according to EC"""
        conc_data = self._load_concrete(file_path='concrete_ec.csv')
        conc_class_names = [el[0] for el in conc_data]
        class_ind = conc_class_names.index(self.cl_conc)
        self._E_cm = float(conc_data[class_ind][9])
        return self._E_cm  # kluczowe pytanie - czym się różni atrybut od właściwości klasy (odp: wła. ma setter i deleter)
    
    def _reinf_heights(self):
        layer_one = self.c + self.fi_s + self.fi * 0.5
        reinf_heights = (0.001 * layer_one, 
                         0.001 * (layer_one + 2 * self.fi),
                         0.001 * (layer_one + 4 * self.fi))
        return reinf_heights
    
    def _initial_rotation(self):
        # ROTACJA TAKA ŻEBY DAŁA  JEDNĄ SETNĄ PROMILA W SKARJNYCH WŁÓKNACH
        fi_init = GeneralAxBend.EPS_INIT / max(self.e_vert, self.h_top)
        return fi_init
    
    def _strains_in_steel(self, eps_cur=0.0, fi_cur=0.0):
        """finds strain in steel layers basing on rotation and eps_linear"""
        r_heights = self._reinf_heights()
        strain_steel = [0] * 6
        # if self.nl_reinf_top == 3:
        #     strain_steel[0] = ((self.h_top - r_heights[0]) * fi_cur)
        #     strain_steel[1] = ((self.h_top - r_heights[1]) * fi_cur)
        for i in range(self.nl_reinf_top):
            strain_steel[i] = - ((self.h_top - r_heights[i]) * fi_cur)
        for i in range(3, 3 + self.nl_reinf_bottom, 1):
            strain_steel[i] = ((self.h_top - r_heights[-i+2]) * fi_cur)
        return strain_steel

    def _stress_strain_conc(self, strain_conc=0.0):
        """returns stress-strain relationsip in cocnrete"""
        f_cd = self._get_fcd_eta()
        if strain_conc >= GeneralAxBend.EC2 and strain_conc <= GeneralAxBend.ECU2:
            sigma_conc = f_cd
        elif strain_conc < GeneralAxBend.EC2 and strain_conc >= 0:
            sigma_conc = f_cd * (1 - (1 - strain_conc / GeneralAxBend.EC2) ** 2)
        elif strain_conc < 0:
            sigma_conc = 0
        else:  # strain_conc > GeneralAxBend.ECU2
            sigma_conc = strain_conc * (f_cd / GeneralAxBend.ECU2)
        return sigma_conc

    def _stress_strain_steel(self, strain_steel=0.0):
        """returns stress-strain relationsip in steel"""
        #k = self.cl_steel_data[6]
        e_uk = self.cl_steel_data[5]
        e_ud = e_uk / 1.15
        f_yd = self.cl_steel_data[1]
        #strain_s_dash = strain_steel / GeneralAxBend.EYD
        #sigma_s_dash = k * strain_s_dash + ((1 - k) * strain_s_dash) / ((1 + strain_s_dash ** GeneralAxBend._R) ** (1 / GeneralAxBend._R))
        #sigma_steel = f_yd * sigma_s_dash
        #k_si = k  + (1 - k) / ((1 + strain_s_dash ** GeneralAxBend._R) ** (1 / GeneralAxBend._R))
        if abs(strain_steel) < GeneralAxBend.EYD:
            sigma_steel = GeneralAxBend.E_STELL * strain_steel
        elif abs(strain_steel) >= GeneralAxBend.EYD and abs(strain_steel) <= e_ud:
            sigma_steel = f_yd
        else:  # abs(strain_steel) > e_yd
            sigma_steel = strain_steel * (f_yd / e_ud)
        # WARINING: doubling of stresses in the reinforcement surrounded 
        # by concrete in compression was not taken into account!!
        return sigma_steel, e_ud
    


    def _find_n_eps_m_eps_fi(self, eps_cur=0.0, fi_cur=0.0):
        """finds corresponding quantities of M N with regard to input
        values of eps_cur, fi_cur
        BY COMPUTING INTEGRALS OF STRESSESS"""
        return m_eps, n_eps, m_fi, e_fi

    def _solve_two_eq(self, eps_cur=0, fi_cur=0):
        """returns new value of eps and fi cur
        by solving linear system of two equations with two variables"""
        # 4x + 3y + 2z = 25
        # -2x + 2y + 3z = -10
        # 3x -5y + 2z = -4
        # A = np.array([[4, 3, 2], [-2, 2, 3], [3, -5, 2]])
        # B = np.array([25, -10, -4])
        # X = np.linalg.inv(A).dot(B)
        # print(X)
        pass
        return eps_new, fi_new

    def optimise_eps_and_fi(self):
        """search iteratively for best values
        of eps and fi so that M_final, and N_fonal were
        close to M_sd and N_sd"""
        return eps_final, fi_final, n_sd_final, m_sd_final

def main():
    my_rc_cross_sec = GeneralAxBend(name='GENERAL_CROSS-SECT_no1',
                                b=(3, 1.5, 4), # [m] width of the individual rectangles
                                h=(1.5, 2.5, 2.5), # [m] height of the individual rectangles
                                hsl=0.20, #[m] thickness of upper slab
                                beff=1.2, #[m] effective width of upper slab
                                cl_conc='C30_37',
                                cl_steel='b500sp',
                                c=30, # [mm]
                                fi=25, # [mm]
                                fi_s=12, # [mm]
                                fi_opp=25, # [mm]
                                nl_reinf_top=3, # [mm] denotes number of layers of top reinforcement
                                nl_reinf_bottom=3, # [mm] denotes number of layers of bottom reinforcement
                                m_sd=11000, # [kNm]
                                n_sd=2000) # [kN]
    print(my_rc_cross_sec._stress_strain_conc(strain_conc=0.0040))
    print(my_rc_cross_sec._stress_strain_steel(0.07))
    print(my_rc_cross_sec.e_vert)
    print(my_rc_cross_sec._reinf_heights())
    print(my_rc_cross_sec.h_top)
    print(my_rc_cross_sec.E_cm)
    print(my_rc_cross_sec._initial_rotation())
    print(my_rc_cross_sec._strains_in_steel(eps_cur=0.0, fi_cur=GeneralAxBend.EPS_INIT))
if __name__ == '__main__':
    main()
