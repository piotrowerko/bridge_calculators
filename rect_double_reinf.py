import csv
import math

from mat_rein_steel import ReinSteel
from rect_single_reinf import RectCrSectSingle


class RectCrSectDoubleR(RectCrSectSingle):
    def __init__(self, name, b, h, cl_conc, cl_steel, c, fi, fi_s, fi_opp):
        super().__init__(name, b, h, cl_conc, cl_steel, c, fi, fi_s)
        self.fi_opp = fi_opp  # diameter of second package of reinforcement
        
    def _compute_a2(self):
        """returns a2 value assuming one row of reinf"""
        return (self.c + 0.5 * self.fi_opp + self.fi_s) / 1000
    
    def _compute_a_s2(self):
        """return an area of e.g. upper reinforcement"""
        num_of_reb = float(input('input num of upper rebars: ').strip())  #.split()
        return math.pi * (self.fi_opp / 1000) ** 2 / 4 * num_of_reb
    
    def _compute_ksi_eff_double_r(self):
        a_s1 = self._compute_a_s1()
        a_s2 = self._compute_a_s2()
        f_yd = self.cl_steel_data[1]
        f_cd = self._get_fcd()
        nominator = a_s1 * f_yd - a_s2 * f_yd
        denominator = self.b * self._compute_d() * f_cd
        return nominator / denominator, a_s1, a_s2
    
    def compute_m_rd_double_r(self):
        ksi_eff, a_s1, a_s2 = self._compute_ksi_eff_double_r()
        d = self._compute_d()
        f_cd = self._get_fcd()
        a2 = self._compute_a2()
        if ksi_eff <= self.cl_steel_data[3]:
            print("reinforcement is fully used; sigma_s = f_yd")
            if ksi_eff > 2 * a2 / d:
                print("ksi_eff > 2 * a2 / d")
                m_rd = 1000 * ksi_eff * (1 - 0.5 * ksi_eff) * d ** 2 * self.b * f_cd \
                + 1000 * (a_s2 * (d - a2) * self.cl_steel_data[1])
                return m_rd, ksi_eff
            else:
                print("ksi_eff <= 2 * a2 / d")
                m_rd = 1000 * a_s1 * (d - a2) * self.cl_steel_data[1]
                return m_rd, ksi_eff
        else:
            print("reinforcement is NOT fully used; sigma_s < f_yd")
            ksi_eff = self.cl_steel_data[3]  # = ksi_eff_lim
            m_rd = 1000 * ksi_eff * (1 - 0.5 * ksi_eff) * d ** 2 * self.b * f_cd
            return m_rd, ksi_eff
    
def main():
    my_double_reinf_cross_sec = RectCrSectDoubleR(name='moj_przekr_prost',
                                                  b=0.5,
                                                  h=1.5,
                                                  cl_conc='C30_37',
                                                  cl_steel='bst500s',
                                                  c=30,
                                                  fi=25,
                                                  fi_s=12,
                                                  fi_opp=25)
    results = my_double_reinf_cross_sec.compute_m_rd_double_r()
    print(f'ksi eff:  {results[1]}')
    print(f'max. load capasity bending moment:  {results[0]}')

if __name__ == '__main__':
    main()
    