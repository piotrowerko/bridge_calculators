import csv
import math

from mat_rein_steel import ReinSteel

class PrzekrProst():
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
    
    def _get_fcd(self, file_path='concrete_ec.csv'):
        conc_data = self._load_concrete(file_path='concerete_ec.csv')
        conc_class_names = [el[0] for el in conc_data]
        class_ind = conc_class_names.index(self.cl_conc)
        return float(conc_data[class_ind][3])
    
    def _compute_a1(self):
        """returns a1 value assuming one row of reinf"""
        return self.c + self.fi + self.fi_s
    
    def _compute_d(self):
        a1 = self._compute_a1()
        return self.h - a1 / 1000
    
    def compute_a_s1(self):
        """return an area of bottom reinforcement"""
        num_of_reb = float(input('input num of rebars: ').strip())  #.split()
        return math.pi * (self.fi / 1000) ** 2 / 4 * num_of_reb
    
    def _compute_ksi_eff(self):
        a_s1 = self.compute_a_s1()
        f_yd = self.cl_steel_data[1]
        f_cd = self._get_fcd()
        nominator = a_s1 * f_yd
        denominator = self.b * self._compute_d() * f_cd
        return nominator / denominator
    
    def _compute_m_rd(self):
        ksi_eff = self._compute_ksi_eff()
        d = self._compute_d()
        f_cd = self._get_fcd()
        if ksi_eff <= self.cl_steel_data[3]:
            print("reinforcement is fully used; sigma_s = f_yd")
            m_rd = 1000 * ksi_eff * (1 - 0.5 * ksi_eff) * d ** 2 * self.b * f_cd
            return m_rd, ksi_eff
        else:
            print("reinforcement is NOT fully used; sigma_s < f_yd")
            ksi_eff = self.cl_steel_data[3]  # = ksi_eff_lim
            m_rd = 1000 * ksi_eff * (1 - 0.5 * ksi_eff) * d ** 2 * self.b * f_cd
            return m_rd, ksi_eff
        
    
def main():
    moj_przekr_prost = PrzekrProst(name='moj_przekr_prost',
                                   b=1,
                                   h=1.5,
                                   cl_conc='C30_37',
                                   cl_steel='bst500s',
                                   c=30,
                                   fi=25,
                                   fi_s=12)
    results = moj_przekr_prost._compute_m_rd()
    print(f'ksi eff:  {results[1]}')
    print(f'max. load capasity bending moment:  {results[0]}')

if __name__ == '__main__':
    main()