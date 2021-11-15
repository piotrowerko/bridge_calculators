class CharGeom():
    def find_center_m(self, a):
        self.a = a
        A = []
        for i in range(0, len(a), 2):
            A.append(a[i] * a[i+1])
        self.A = A
        ex = (A[0] * 0.5 * a[1] + A[1] * (a[1] + 0.5 * a[3]) 
              + A[2] * (a[1] + a[3] + 0.5 * a[5])) / (sum(A))
        self.ex = ex
        return ex


def main():
    moje_cw = CharGeom()
    moje_cw.find_center_m((3, 2, 1.5, 3, 5, 2))
    print(f'pole: {moje_cw.A}; środek ciężkości: {moje_cw.ex}')


if __name__ == '__main__':
    main()