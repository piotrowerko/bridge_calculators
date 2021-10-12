class PrzekrKabl():
    def __init__(self, name, rozstaw, l1, l2, bdz, hpl):  
        """l1 - długośćprzęsła skrajnego, 
        l2 - centralnego, 
        bdz - szerokośćsrednia dźwiagara
        hpl - grubosc płyty pomostowej"""
        self.name = name
        self.rozstaw = rozstaw
        self.l1 = l1
        self.l2 = l2
        self.bdz = bdz
        self.hpl = hpl
    
   
    def __str__(self):
        return f'{self.name}'
    
    def compute_beff(self):
        """wyznacza pełną szerokość współpracującą dźwigara betonowego"""
        return (0.2 * self.l2 + self.bdz, self.hpl * 12 + self.bdz, self.rozstaw)

def main():
    moj_przekr_trap = PrzekrKabl('moj_przekr_trap', 4, 0, 30, 1.4, 0.25)
    pp = moj_przekr_trap.compute_beff()
    print(f'minimum z: {pp}')
    print(f'to: beff = {min(pp)}')
    

if __name__ == '__main__':
    main()