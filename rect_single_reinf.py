import csv
import math

from mat_rein_steel import ReinSteel

class RectCrSectSingle():
    """bending of rectangular cross section with singular reinforcement:
    evaluation of bending moment capasity [kNm]"""
    def __init__(self, name, b, h, cl_conc, cl_steel, c, fi, fi_s):  
        self.name = name
        self.b = b  # cross section width
        self.h = h  # cross section height
        self.cl_conc = cl_conc
        self.cl_steel_data = ReinSteel.REIN_STEEL_DATA[cl_steel]
        self.c = c
        self.fi = fi
        self.fi_s = fi_s  # stirrup diameter 

    def __str__(self):
        return f'{self.name}'

    def _load_concrete(self, file_path='concrete_ec.csv'):
        with open(file_path) as fp:
            reader = csv.reader(fp, delimiter=";", quotechar='"')
            _data = [row for row in reader]
            data = list(map(list, zip(*_data)))
        return data
    
    def _get_fcd_eta(self):
        """returns f_cd * eta value according to EC"""
        conc_data = self._load_concrete(file_path='concrete_ec.csv')
        conc_class_names = [el[0] for el in conc_data]
        class_ind = conc_class_names.index(self.cl_conc)
        fck = float(conc_data[class_ind][1])
        fcd = float(conc_data[class_ind][3])
        if fck <= 50:
            eta = 1.0  # effective compressive strength coefficient
        else:
            eta = 1.0 - (fck - 50) / 200
        return eta * fcd
    
    def _compute_a1(self):
        """returns a1 value assuming one row of reinf"""
        return (self.c + 0.5 * self.fi + self.fi_s) / 1000
    
    def _compute_d(self):
        a1 = self._compute_a1()
        return self.h - a1
    
    def _compute_a_s1(self):
        """return an area of e.g. bottom reinforcement"""
        num_of_reb = float(input('input num of bootom rebars: ').strip())  #.split()
        return math.pi * (self.fi / 1000) ** 2 / 4 * num_of_reb
    
    def _compute_ksi_eff_single_r(self):
        a_s1 = self._compute_a_s1()
        f_yd = self.cl_steel_data[1]
        f_cd = self._get_fcd_eta()
        nominator = a_s1 * f_yd
        denominator = self.b * self._compute_d() * f_cd
        return nominator / denominator
    
    def compute_m_rd_single_r(self):
        ksi_eff = self._compute_ksi_eff_single_r()
        d = self._compute_d()
        f_cd = self._get_fcd_eta()
        if ksi_eff <= self.cl_steel_data[3]:  # ksi_eff_lim assuming LAMBDA = 0.8 and EPSILON CU = 3.5 * (10 ** -3)
            print("reinforcement is fully used; sigma_s = f_yd")
            m_rd = 1000 * ksi_eff * (1 - 0.5 * ksi_eff) * d ** 2 * self.b * f_cd
            return m_rd, ksi_eff
        else:
            print("reinforcement is NOT fully used; sigma_s < f_yd")
            ksi_eff = self.cl_steel_data[3]  # = ksi_eff_lim
            m_rd = 1000 * ksi_eff * (1 - 0.5 * ksi_eff) * d ** 2 * self.b * f_cd
            return m_rd, ksi_eff
        
    
def main():
    moj_przekr_prost = RectCrSectSingle(name='moj_przekr_prost',
                                   b=0.5,
                                   h=1.5,
                                   cl_conc='C30_37',
                                   cl_steel='bst500s',
                                   c=30,
                                   fi=25,
                                   fi_s=12)
    results = moj_przekr_prost.compute_m_rd_single_r()
    print(f'ksi eff:  {results[1]}')
    print(f'max. load capasity bending moment:  {results[0]}')

if __name__ == '__main__':
    main()