class DwaWektory():
    def __init__(self, Ax, Ay, Bx, By, Cx, Cy):
        self.AB = [Bx - Ax, By - Ay]
        self.AC = [Cx - Ax, Cy - Ay]
        self.lAB = (self.AB[0] ** 2 + self.AB[1] ** 2) ** 0.5
        self.lAC = (self.AC[0] ** 2 + self.AC[1] ** 2) ** 0.5
        self.ABoAC = self.AB[0] * self.AC[0] + self.AB[1] * self.AC[1]


def main():
    moja_para_wektorow = DwaWektory(3, 7, 3, 3, 7, 3)
    print(f'wsp AB = {moja_para_wektorow.AB}')
    print(f'wsp AC = {moja_para_wektorow.AC}')
    print(f'długość AB = {moja_para_wektorow.lAB}')
    print(f'długość AC = {moja_para_wektorow.lAC}')
    print(f'iloczyn skalarny = {moja_para_wektorow.ABoAC}')
    
    moja_para_wektorow2 = DwaWektory(2, 8, 2, 2, 8, 2)
    print(f'wsp AB = {moja_para_wektorow2.AB}')
    print(f'wsp AC = {moja_para_wektorow2.AC}')
    print(f'długość AB = {moja_para_wektorow2.lAB}')
    print(f'długość AC = {moja_para_wektorow2.lAC}')
    print(f'iloczyn skalarny = {moja_para_wektorow2.ABoAC}')
    

if __name__ == '__main__':
    main()