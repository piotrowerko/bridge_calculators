"""
sources of knoweladge
pietryga
https://chodor-projekt.net/encyclopedia/krzywe-interakcji-m-n-zelbetu/
knauf
łapko jensen
oleszek
"""

import csv
import math
import numpy as np

from t_sect_ben_reinf import TCrReinf

class GeneralAxBend(TCrReinf):
    """cross section axial force + bending moment capasity computations;
    returns stresses in conrete and steel layers according to computed itaratively
    (gradient descent) epsilon and fi (strain and rotation of cross section) 
    to achive close resulting M and N as M_sd and N_sd"""
    EC2 = 0.0020  # concrete strain of section in the begining of flat (plastic) stresses
    ECU2 = 0.0035  # max. concrete strain of section
    EYD = 0.00217  # max. steel strain (FOR BOTH B500SP AND BST500S STEEL)
    _R = 20  # max. steel strain (FOR BOTH B500SP AND BST500S STEEL ?)
    
    def __init__(self, name, b, h, hsl, beff, cl_conc, cl_steel, c, fi, fi_s, fi_opp, m_sd, n_sd):
        super().__init__(name, b, h, hsl, beff, cl_conc, cl_steel, c, fi, fi_s, fi_opp, m_sd)
        self.n_sd = n_sd
    
    def _get_E_cm(self):
        """returns E_cm value according to EC"""
        conc_data = self._load_concrete(file_path='concrete_ec.csv')
        conc_class_names = [el[0] for el in conc_data]
        class_ind = conc_class_names.index(self.cl_conc)
        E_cm = float(conc_data[class_ind][9])
        return E_cm
    
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
            sigma_conc = 0
        return sigma_conc

    def _stress_strain_steel(self, strain_steel=0.0):
        """returns stress-strain relationsip in steel"""
        k = self.cl_steel_data[6]
        f_yd = self.cl_steel_data[1]
        strain_s_dash = strain_steel / GeneralAxBend.EYD
        #sigma_s_dash = k * strain_s_dash + ((1 - k) * strain_s_dash) / ((1 + strain_s_dash ** GeneralAxBend._R) ** (1 / GeneralAxBend._R))
        #sigma_steel = f_yd * sigma_s_dash
        k_si = k  + (1 - k) / ((1 + strain_s_dash ** GeneralAxBend._R) ** (1 / GeneralAxBend._R))
        return k_si

    def _find_n_eps_m_eps_fi(self, eps_cur=0, fi_cur=0):
        """finds corresponding quantities of M N with regard to input
        values of eps_cur, fi_cur
        BY COMPUTING INTEGRALS OF STRESSESS"""
        return m_eps, n_eps, m_fi, e_fi

    def _solve_two_eq(self, eps_cur=0, fi_cur=0):
        """returns new value of eps and fi cur
        by solving linear system of two equations with two wariables"""
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
                                b=0.5, # [m]
                                h=1.5, # [m]
                                hsl=0.20, #[m] thickness of upper slab
                                beff=1.2, #[m] effective width of upper slab
                                cl_conc='C30_37',
                                cl_steel='b500sp',
                                c=30, # [mm]
                                fi=25, # [mm]
                                fi_s=12, # [mm]
                                fi_opp=25, # [mm]
                                m_sd=11000, # [kNm]
                                n_sd=2000) # [kN]
    res = my_rc_cross_sec._get_E_cm()
    print(res)
    print(my_rc_cross_sec._stress_strain_conc(strain_conc=0.0000))
    print(my_rc_cross_sec._stress_strain_steel(0.05))


if __name__ == '__main__':
    main()