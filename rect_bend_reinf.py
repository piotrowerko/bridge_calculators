import csv
import math

from mat_rein_steel import ReinSteel
from rect_single_reinf import RectCrSectSingle
from rect_double_reinf import RectCrSectDoubleR


class RectCrReinf(RectCrSectDoubleR):
    def __init__(self, name, b, h, cl_conc, cl_steel, c, fi, fi_s, fi_opp, m_sd):
        super().__init__(name, b, h, cl_conc, cl_steel, c, fi, fi_s, fi_opp)
        self.m_sd = m_sd

    def __str__(self):
        return f'{self.name}'

    @staticmethod
    def _equationroots(_b, c, a=1):
        """computes roots of quadratic equation""" 
      # calculating discriminant using formula
        dis = _b * _b - 4 * a * c 
        sqrt_val = math.sqrt(abs(dis)) 
        # checking condition for discriminant
        if dis > 0:  # real and different roots 
            r1 = (-_b + sqrt_val)/(2 * a)
            r2 = (-_b - sqrt_val)/(2 * a)
        elif dis == 0:  # real and same roots
            r1 = -_b / (2 * a)
            r2 = None
        # when discriminant is less than 0
        else:  # Complex Roots
            r1 = (- _b / (2 * a), " + i", sqrt_val) 
            r2 = (- _b / (2 * a), " - i", sqrt_val)
        return r1, r2

    def _compute_ksi_eff_reinf(self, b=None):
        if b == None:
            b = self.b
        f_cd = self._get_fcd()
        d = self._compute_d()
        c = 2 * (self.m_sd * 0.001) / (b * d ** 2 * f_cd)  # "c" here is NOT concrete cover
        return RectCrReinf._equationroots(_b=-2, c=c, a=1)

    def compute_reinf_rect(self, b=None):
        if b == None:
            b = self.b
        d = self._compute_d()
        f_cd = self._get_fcd()
        f_yd = self.cl_steel_data[1]
        ksi_eff_lim = self.cl_steel_data[3]
        a2 = self._compute_a2()
        ksi_eff_two = self._compute_ksi_eff_reinf(b)
        print(f'ksi eff roots: {ksi_eff_two}')
        ksi_eff = min(ksi_eff_two)
        print(ksi_eff)
        if ksi_eff > ksi_eff_lim:
            print('As2 reinforcement is necessary')
            ksi_eff = ksi_eff_lim  # = ksi_eff_lim
            m_rd_star = 1000 * ksi_eff * (1 - 0.5 * ksi_eff) * d ** 2 * b * f_cd
            delta_m = self.m_sd - m_rd_star
            a_s1_star = ksi_eff_lim * d * b * f_cd / f_yd
            a_s2 = (delta_m / 1000) / (f_yd * (d - a2))  # = As1 star star
            a_s1 = a_s1_star + a_s2
        else:
            print('only the As1 reinforcement is sufficient')
            a_s1 = ksi_eff * d * b * f_cd / f_yd
            a_s2 = 0
        n_s1 = math.ceil(a_s1 / (math.pi * (self.fi / 1000) ** 2 / 4 ))
        n_s2 = math.ceil(a_s2 / (math.pi * (self.fi / 1000) ** 2 / 4 ))
        return a_s1, n_s1, a_s2, n_s2
    
def main():
    my_rc_cross_sec = RectCrReinf(name='moj_przekr_prost',
                                b=1.2, # [m]
                                h=1.5, # [m]
                                cl_conc='C30_37',
                                cl_steel='bst500s',
                                c=30, # [mm]
                                fi=25, # [mm]
                                fi_s=12, # [mm]
                                fi_opp=25, # [mm]
                                m_sd=5000) # [kNm]
    res = my_rc_cross_sec.compute_reinf_rect()
    print(f'bottom needs: {res[1]} x fi{my_rc_cross_sec.fi} >= {res[0]}')
    print(f'upper needs: {res[3]} x fi{my_rc_cross_sec.fi} >= {res[2]}')

if __name__ == '__main__':
    main()