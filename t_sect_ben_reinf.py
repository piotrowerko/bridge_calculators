import csv
import math

from mat_rein_steel import ReinSteel
from rect_single_reinf import RectCrSectSingle
from rect_double_reinf import RectCrSectDoubleR
from rect_bend_reinf import RectCrReinf

class TCrReinf(RectCrReinf):
    """T cross section bending computations"""
    
    def __init__(self, name, b, h, hsl, beff, cl_conc, cl_steel, c, fi, fi_s, fi_opp, m_sd):
        super().__init__(name, b, h, cl_conc, cl_steel, c, fi, fi_s, fi_opp, m_sd)
        self.hsl = hsl
        self.beff = beff

    def __str__(self):
        return f'{self.name}'

    def _compute_beta(self):
        d = self._compute_d()
        return self.hsl / d, d

    def _compute_mrd_sl(self):
        """returns load capasity of upper slab Mrd [kNm]"""
        beta, d = self._compute_beta()
        f_cd = self._get_fcd()
        m_rd = beta * (1 - 0.5 * beta) * d ** 2 * self.beff * f_cd * 1000
        return m_rd, beta, d, f_cd

    def _compute_ksi_eff_reinf_T(self, m_rd_1):
        # a_s1 = self._compute_a_s1()
        # f_yd = self.cl_steel_data[1]
        f_cd = self._get_fcd()
        d = self._compute_d()
        c = 2 * 0.001 * (self.m_sd - m_rd_1) / (self.b * d ** 2 * f_cd)
        return RectCrReinf._equationroots(b=-2, c=c, a=1)

    def compute_reinf(self):
        m_rd, beta, d, f_cd = self._compute_mrd_sl()
        f_yd = self.cl_steel_data[1]
        ksi_eff_lim = self.cl_steel_data[3]
        if m_rd < self.m_sd:
            a_s11 = beta * d * (self.beff - self.b) * f_cd / f_yd
            # computes load capasity of uuper slab without web [kNm]:
            m_rd_1 = beta * (1 - 0.5 * beta) * d ** 2 * (self.beff - self.b) * f_cd
            ksi_eff_two = self._compute_ksi_eff_reinf_T(m_rd_1)
            print(f'ksi eff roots: {ksi_eff_two}')
            ksi_eff = min(ksi_eff_two)
            if ksi_eff > ksi_eff_lim:
                pass

def main():
    my_rc_cross_sec = TCrReinf(name='moj_przekr_T',
                                b=0.5, # [m]
                                h=1.5, # [m]
                                hsl=0.20, #[m] thickness of upper slab
                                beff=1.2, #[m] effective width od upper slab
                                cl_conc='C30_37',
                                cl_steel='bst500s',
                                c=30, # [mm]
                                fi=25, # [mm]
                                fi_s=12, # [mm]
                                fi_opp=25, # [mm]
                                m_sd=11000) # [kNm]
    print(f'{my_rc_cross_sec.compute_reinf()}')

if __name__ == '__main__':
    main()
