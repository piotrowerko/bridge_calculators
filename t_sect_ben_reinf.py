import csv
import math

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
        m_rd = 1000 * beta * (1 - 0.5 * beta) * d ** 2 * self.beff * f_cd
        return m_rd, beta, d, f_cd

    def _compute_ksi_eff_reinf_T(self, m_rd_1):
        # a_s1 = self._compute_a_s1()
        # f_yd = self.cl_steel_data[1]
        f_cd = self._get_fcd()
        d = self._compute_d()
        c = 2 * 0.001 * (self.m_sd - m_rd_1) / (self.b * d ** 2 * f_cd)
        return RectCrReinf._equationroots(_b=-2, c=c, a=1)

    def compute_reinf_T(self):
        m_rd, beta, d, f_cd = self._compute_mrd_sl()
        f_yd = self.cl_steel_data[1]
        ksi_eff_lim = self.cl_steel_data[3]
        a2 = self._compute_a2()
        if m_rd < self.m_sd:  # the compression zone extends below the plate
            print('the compression zone extends below the plate')
            a_s11 = beta * d * (self.beff - self.b) * f_cd / f_yd
            # computes load capasity of uuper slab without web [kNm]:
            m_rd_1 = 1000 * beta * (1 - 0.5 * beta) * d ** 2 * (self.beff - self.b) * f_cd
            print(f'Mrd1: {m_rd_1}')
            ksi_eff_two = self._compute_ksi_eff_reinf_T(m_rd_1)
            print(f'ksi eff roots: {ksi_eff_two}')
            ksi_eff = min(ksi_eff_two)
            if ksi_eff <= ksi_eff_lim:
                print(f'only the As1 reinforcement is sufficient {ksi_eff}, {ksi_eff_lim}')
                a_s12 = ksi_eff * d * self.b * f_cd / f_yd
                a_s2 = 0
                a_s1 = a_s11 + a_s12
            else:
                print('As2 reinforcement is necessary')
                ksi_eff = ksi_eff_lim  # = ksi_eff_lim
                a_s12 = ksi_eff * d * self.b * f_cd / f_yd
                # computation of limit load capasity of the web:
                m_rd_web = 1000 * ksi_eff * (1 - 0.5 * ksi_eff) * d ** 2 * self.b * f_cd
                # computation of additional compression and tension reinforcement:
                a_s13 = 0.001 * (self.m_sd - m_rd_1 - m_rd_web) / ((d - a2) * f_yd)
                a_s2 = a_s13
                a_s1 = a_s11 + a_s12 + a_s13
            n_s1 = math.ceil(a_s1 / (math.pi * (self.fi / 1000) ** 2 / 4 ))
            n_s2 = math.ceil(a_s2 / (math.pi * (self.fi / 1000) ** 2 / 4 ))
            return a_s1, n_s1, a_s2, n_s2
        else: # the compression zone does NOT extends below the plate
            print('the compression zone does NOT extends below the plate')
            return self.compute_reinf_rect(b=self.beff)

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
    res = my_rc_cross_sec.compute_reinf_T()
    print(f'bottom needs: {res[1]} x fi{my_rc_cross_sec.fi} >= {res[0]}')
    print(f'upper needs: {res[3]} x fi{my_rc_cross_sec.fi} >= {res[2]}')

if __name__ == '__main__':
    main()
