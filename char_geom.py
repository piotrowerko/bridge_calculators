class CharGeom():
    def find_center_m(self, a):
        areas = []
        for i in range(0, len(a), 2):
            areas.append(a[i] * a[i+1])
        ex = (areas[0] * 0.5 * a[1] + areas[1] * (a[1] + 0.5 * a[3]) 
              + areas[2] * (a[1] + a[3] + 0.5 * a[5])) / (sum(areas))
        return ex, areas, sum(areas) 

    def mom_of_int_sin(self, a):
        mom_of_int_sin = []
        for i in range(0, len(a), 2):
            mom_of_int_sin.append(a[i] * a[i+1] ** 3 / 12)
        return mom_of_int_sin

    def mom_of_int_bend(self, a):
        ex, areas, area = self.find_center_m(a)
        mom_of_int_sin = self.mom_of_int_sin(a)
        mom_of_int_bend = sum(mom_of_int_sin) +  areas[0] * abs((ex - 0.5 * a[1])) ** 2 \
            + areas[1] * abs((ex - (a[1] + 0.5 * a[3]))) ** 2 \
            + areas[2] * abs((a[1] + a[3] + 0.5 * a[5] - ex)) ** 2
        return mom_of_int_bend

def main():
    moje_cw = CharGeom()
    a = (10, 0.1, 0.70, 0.3, 0, 0)
    res = moje_cw.find_center_m(a)
    print(f'pole: {res[1]}; pole razem: {res[2]}; środek ciężkości: {res[0]}')
    #print(sum(moje_cw.mom_of_int_sin(a)))
    print(f'mom bezw: {moje_cw.mom_of_int_bend(a)}')

if __name__ == '__main__':
    main()