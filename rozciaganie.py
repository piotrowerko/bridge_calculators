class PretRozc():
    def __init__(self, l, a, b, E, P):
        self.A = a * b / (100 * 100)
        self.sigma = (P / 1000) / self.A
        self.dl = self.sigma * l * 1000 / (E * 1000)
        

def main():
    pret1 = PretRozc(4.80, 2, 7, 190, 380)
    pret2 = PretRozc(4.00, 3, 5, 215, 310)
    print(f'naprężenie w [MPa]: {pret1.sigma}')
    print(f'wydłużenie w [mm]: {pret1.dl}')
    print(f'naprężenie 2 w [MPa]: {pret2.sigma}')
    print(f'wydłużenie 2 w [mm]: {pret2.dl}')


if __name__ == '__main__':
    main()