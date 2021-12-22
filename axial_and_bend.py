"""
sources of knoweladge
knauf
łapko jensen
pietryga
https://chodor-projekt.net/encyclopedia/krzywe-interakcji-m-n-zelbetu/
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
    EPS_INIT = 0.001 * 10 ** -3
    N_CONC_LAYERS = 10  # number of layers of virtual division of the concrete cross-section for the needs of numerical integrals
    CONC_L_THIC = 0.1 # concrete layer thickness of virtual division of the concrete cross-section for the needs of numerical integrals
    

    def __init__(self, name, b, h, hsl, 
                 beff, cl_conc, cl_steel, 
                 c, fi, fi_s, fi_opp, 
                 nl_reinf_top, nl_reinf_bottom,
                 m_sd, n_sd):
        super().__init__(name, b, h, hsl, beff, cl_conc, cl_steel, c, fi, fi_s, fi_opp, m_sd)
        self.n_sd = n_sd
        char_geom = CharGeom()
        self.e_vert = char_geom.find_center_m((b[0], h[0], b[1], h[1], b[2], h[2]))[0]  # = sefl.h_bottom
        self.h_top = sum(self.h) - self.e_vert
        self.nl_reinf_top = nl_reinf_top
        self.nl_reinf_bottom = nl_reinf_bottom
        
    @classmethod
    def alter_constr(self):
        """alternative contructor with custom reinforcement heights"""
        pass
    
    @property
    def E_cm(self):
        """returns E_cm value according to EC"""
        conc_data = self._load_concrete(file_path='concrete_ec.csv')
        conc_class_names = [el[0] for el in conc_data]
        class_ind = conc_class_names.index(self.cl_conc)
        self._E_cm = float(conc_data[class_ind][9])
        return self._E_cm  # kluczowe pytanie - czym się różni atrybut od właściwości klasy (odp: wła. ma setter i deleter)

    def _initial_rotation(self):
        # ROTACJA TAKA ŻEBY DAŁA  JEDNĄ SETNĄ PROMILA W SKRAJNYCH WŁÓKNACH
        fi_init = GeneralAxBend.EPS_INIT / max(self.e_vert, self.h_top)
        return fi_init
    
    def _conc_layers(self):
        """finds heights and number of layers of virtual division in the concrete cross-section"""
        n_layers = int(sum(self.h) / GeneralAxBend.CONC_L_THIC)
        conc_lay_heights = []
        conc_lay_widths = []
        conc_lay_areas = []
        for i in range(n_layers):
            curr_rel_height = self.h_top - (0.5 + i) * GeneralAxBend.CONC_L_THIC
            curr_height = sum(self.h) - (0.5 + i) * GeneralAxBend.CONC_L_THIC
            conc_lay_heights.append(curr_rel_height)
            if curr_height < self.h[0]:
                conc_lay_widths.append(self.b[0])
            elif curr_height > self.h[0] and curr_height < (self.h[0]+self.h[1]):
                conc_lay_widths.append(self.b[1])
            else:
                conc_lay_widths.append(self.b[2])
            conc_lay_areas.append(conc_lay_heights[i] * conc_lay_widths[i])
        return n_layers, conc_lay_heights, conc_lay_widths, conc_lay_areas
    
    def _strains_in_conc(self, eps_cur=0.0, fi_cur=0.0):
        """finds strain in concrete layers basing on rotation and eps_linear"""
        n_layers, conc_lay_heights, conc_lay_widths, conc_lay_areas = self._conc_layers()
        strain_conc = [(conc_lay_heights[i] * fi_cur + eps_cur) for i in range(n_layers)]
        return strain_conc, n_layers, conc_lay_heights, conc_lay_widths, conc_lay_areas
    
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
            
    def _stress_in_conc(self, eps_cur=0.0, fi_cur=0.0):
        """returns stresses in concrete layers"""
        strain_conc, n_layers, conc_lay_heights, conc_lay_widths, conc_lay_areas = self._strains_in_conc(eps_cur, fi_cur)
        stress_conc = [self._stress_strain_conc(el) for el in strain_conc]
        return stress_conc, n_layers, conc_lay_heights, conc_lay_widths, conc_lay_areas
    
    def _reinf_geom(self):
        layer_one = self.c + self.fi_s + self.fi * 0.5
        n_bott, rebar_numbers_bott = self.nl_reinf_bottom
        n_upp, rebar_numbers_upp = self.nl_reinf_top
        n_of_layers = n_bott + n_upp
        reinf_heights = [(0.001 * layer_one + 0.002 * i * self.fi) for i in range(n_bott)] \
            + [(sum(self.h) - 0.001 * layer_one - 0.002 * i * self.fi) for i in range(n_upp)]
        single_bar_area =  math.pi * ((0.001 * self.fi) ** 2) / 4
        reinf_areas = [rebar_numbers_bott[i] * single_bar_area for i in range(n_bott)] \
            + [rebar_numbers_upp[i] * single_bar_area for i in range(n_upp)]
        return reinf_heights, reinf_areas, n_of_layers
    
    def _strains_in_steel(self, eps_cur=0.0, fi_cur=0.0):
        """finds strain in steel layers basing on rotation and eps_linear"""
        r_heights, reinf_areas, n_of_layers = self._reinf_geom()
        strain_steel = [0] * n_of_layers
        for i in range(self.nl_reinf_top[0]):
            strain_steel[i] = ((self.h_top - r_heights[i]) * fi_cur) + eps_cur
        for i in range(3, 3 + self.nl_reinf_bottom[0], 1):
            strain_steel[i] = ((self.e_vert - r_heights[-i+2]) * fi_cur) + eps_cur
        return strain_steel, r_heights, reinf_areas
    
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
    
    def _stress_in_steel(self, eps_cur=0.0, fi_cur=0.0):
        """returns stresses in reinforcement layers"""
        strain_steel, r_heights, reinf_areas = self._strains_in_steel(eps_cur, fi_cur)
        stress_steel = [self._stress_strain_steel(el)[0] for el in strain_steel]
        return stress_steel, r_heights, reinf_areas
    
    def _relative_heights(self, heights):
        "returns heights relative to center of gravity"
        rel_heights = [(i - self.e_vert) for i in heights]
        return rel_heights
    
    def _bending_moment(self):
        """numerical integral for internal bending moment computation;
        takes given strain and rotation,
        returns corresponding bending moment"""
        # suma pole_plastra*naprezenieplastra*odlegloscdoSRC
        # suma pole_warw_pretow*naprezenia_w_wawie-pretow*odlegl_doSRC
        pass
    
    def _axial_force(self, eps_cur, fi_cur):
        """numerical integral for internal axial force computation;
        takes given strain and rotation,
        returns corresponding axial force"""
        stress_in_conc, conc_lay_heights, conc_areas = self._stress_in_conc(eps_cur, fi_cur)[0], \
        self._stress_in_conc(eps_cur, fi_cur)[2], \
        self._stress_in_conc(eps_cur, fi_cur)[4]
        
        stress_in_steel, r_heights, steel_areas = self._stress_in_steel(eps_cur, fi_cur)
        # suma pole_plastra*naprezenieplastra
        # suma pole_warw_pretow*naprezenia_w_wawie-pretow
        forces_in_conc = [stress_in_conc[i] * conc_areas[i] \
            for i in range(len(stress_in_conc))]
        forces_in_steel = [stress_in_steel[i] * steel_areas[i] \
            for i in range(len(stress_in_steel))]
        axial_force = sum(forces_in_conc + forces_in_steel)
        # for i in range(self.nl_reinf_top[0]):
        #     strain_steel[i] = ((self.h_top - r_heights[i]) * fi_cur) + eps_cur
        # for i in range(3, 3 + self.nl_reinf_bottom[0], 1):
        #     strain_steel[i] = ((self.e_vert - r_heights[-i+2]) * fi_cur) + eps_cur
        rel_heights = self._relative_heights(conc_lay_heights)
        bending_moment_conc = [forces_in_conc[i] * rel_heights[i] for i in range(len(forces_in_conc))]
        #DOPISZ SIŁĘ OSIOWĄ !!
        return axial_force, sum(bending_moment_conc)

    def _find_n_eps_m_eps_fi(self, eps_cur=0.0, fi_cur=0.0):
        """finds corresponding quantities of M N with regard to input
        values of eps_cur, fi_cur
        BY COMPUTING INTEGRALS OF STRESSESS"""
        pass
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
        of eps and fi so that M_final, and N_final were
        close to M_sd and N_sd"""
        pass
        return eps_final, fi_final, n_sd_final, m_sd_final

def main():
    my_rc_cross_sec = GeneralAxBend(name='GENERAL_CROSS-SECT_no1',
                                b=(1, 0.3, 1), # [m] width of the individual rectangles
                                h=(0.3, 1, 0.3), # [m] height of the individual rectangles
                                hsl=0.20, #[m] thickness of upper slab
                                beff=1.2, #[m] effective width of upper slab
                                cl_conc='C30_37',
                                cl_steel='b500sp',
                                c=30, # [mm]
                                fi=25, # [mm]
                                fi_s=12, # [mm]
                                fi_opp=25, # [mm]
                                nl_reinf_top=(3, (6, 5, 4)), # [mm] denotes number of layers of top reinforcement and corresponding numbers of rebars
                                nl_reinf_bottom=(3, (6, 5 , 4)), # [mm] denotes number of layers of bottom reinforcement and corresponding numbers of rebars
                                m_sd=11000, # [kNm]
                                n_sd=2000) # [kN]
    # print(my_rc_cross_sec._stress_strain_conc(strain_conc=0.0040))
    # print(my_rc_cross_sec._stress_strain_steel(0.07))
    # print(my_rc_cross_sec.e_vert)
    # print(my_rc_cross_sec._reinf_geom())
    # print(my_rc_cross_sec.h_top)
    # print(my_rc_cross_sec.E_cm)
    init_fi = my_rc_cross_sec._initial_rotation()
    # print(init_fi)
    print('strain in steel', my_rc_cross_sec._strains_in_steel(eps_cur=0.000000, fi_cur=-init_fi)[0])
    # # print(my_rc_cross_sec._conc_layers())
    print('strain in conc',my_rc_cross_sec._strains_in_conc(eps_cur=0.000000, fi_cur=-init_fi)[0])
    print('str in conc:', my_rc_cross_sec._stress_in_conc(eps_cur=0.000000, fi_cur=-init_fi)[0])
    print('str in steel:', my_rc_cross_sec._stress_in_steel(eps_cur=0.000000, fi_cur=-init_fi)[0])
    print('axial force:', my_rc_cross_sec._axial_force(eps_cur=0.000000, fi_cur=-init_fi))
if __name__ == '__main__':
    main()
