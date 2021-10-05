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
    
    def compute_be1_betot(self):
        return (0.7 * self.l2 / 8, (self.rozstaw - self.b0) / 2)

def main():
    moj_przekr_zesp = PrzekrZesp('moj_przekr_zesp', 2.15, 0, 30, 0.15)
    a, b = moj_przekr_zesp.compute_be1_betot()
    print(f'minimum z: {a, b}')
    print(f'to: {min(a, b)}')
    

if __name__ == '__main__':
    main()