class PrzekrZesp():
    def __init__(self, name, rozstaw, l1, l2, b0):  
        """l1 długośćprzęsła skrajnego, 
        l2 centralnego, 
        b0 pomiędzy zewnętrznymi sworzniami zespolajacymi"""
        self.name = name
        self.rozstaw = rozstaw
        self.l1 = l1
        self.l2 = l2
        self.b0 = b0
    
    def __str__(self):
        return f'{self.name}'
    
    def compute_be1_prz(self):
        """wyznacza be1
        a więc odległość efektywną od sworznia do granicy
        współpracy płyty z jednej strony; w przekroju przeslowym"""
        return (0.7 * self.l2 / 8, (self.rozstaw - self.b0) / 2)
    
    def compute_be1_podp(self):
        """wyznacza be1
        a więc odległość efektywną od sworznia do granicy
        współpracy płyty z jednej strony; w podporowym"""
        return (0.25 * (self.l2 + self.l1) / 8, (self.rozstaw - self.b0) / 2)  

def main():
    moj_przekr_zesp = PrzekrZesp('moj_przekr_zesp', 2.15, 20, 30, 0.15)
    pp = moj_przekr_zesp.compute_be1_prz()
    print(f'minimum z: {pp}')
    print(f'to: be1 (przeslo)= {min(pp)}')
    pp1 = moj_przekr_zesp.compute_be1_podp()
    print(f'minimum z: {pp1}')
    print(f'to: be1 (podpora) = {min(pp1)}')
    

if __name__ == '__main__':
    main()