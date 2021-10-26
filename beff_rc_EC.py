class PrzekrKabl():
    def __init__(self, name, rozstaw, li, bdz, hpl, b1, b2):  
        """li - krotka z długościami przęsęł,
        bdz - szerokość dźwiagara w połaczeniu z płytą
        hpl - grubosc płyty pomostowej
        b1,b2, lewa i prawa długość wysięgu półki do sąsiadniego
        dźwigara lub do końca wspornika w przekr. poprzecznym"""
        self.name = name
        self.rozstaw = rozstaw
        self.li = li
        self.bdz = bdz
        self.hpl = hpl
        self.b1 = b1
        self.b2 = b2

    def __str__(self):
        return f'{self.name}'

    def _find_l0(self):
        """wyznacza zestaw długości l0
        źródło: rys. 5.2 EC2-1"""
        l0 = {
            "l0_prz_skr": 0.85 * self.li[0], 
            "l0_podp1" : 0.15 * (self.li[0]+self.li[1]),
            "l0_prz_wewn" : 0.70 * (0.70 * self.li[1]),
            "l0_wsprnik3" : 0.15 * self.li[1] + self.li[2]
        }
        return l0

    def compute_beff_EC(self):
        """wyznacza szerokości współpracującą dźwigara betonowego
        źródło: rys. 5.3 EC2-1"""
        l0 = self._find_l0()
        beffs = []
        keys = []
        for key, val in l0.items():
            keys.append('beff'+key[2:])
            beffs.append([min(0.2 * self.b1 + 0.1 * val, 0.2 * val, self.b1), 
                          min(0.2 * self.b2 + 0.1 * val, 0.2 * val, self.b2)])
        beff12_dict = dict(zip(keys, beffs))
        beffs_tot = []
        for val in beffs:
            beffs_tot.append(val[0]+val[1]+self.bdz)
        beff_tot_dict = dict(zip(keys, beffs_tot))
        return beff12_dict, beff_tot_dict #(0.2 * self.li[1] + self.bdz, self.hpl * 12 + self.bdz, self.rozstaw)

def main():
    moj_przekr_trap = PrzekrKabl('moj_przekr_trap', 4, (20, 30, 20), 1.4, 0.25, 1.8, 1.9)
    pp = moj_przekr_trap.compute_beff_EC()
    print(f'beffs 1, 2 (lewa, prawa): {pp[0]}')
    print(f'beffs total: {pp[1]}')

if __name__ == '__main__':
    main()