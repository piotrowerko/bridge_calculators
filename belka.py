class Belka():
    def __init__(self, a, b, c):
        self.a = a
        self.b = b
        self.c = c
    
    
    def compute_reactions(self, P):
        self.P = P
        self.Rb = (P[0] * (self.a) - P[1] * self.b[0] - 
                   P[2] * sum(self.b[0:2]) - 
                           P[3] * (sum(self.b) + self.c)) / sum(self.b)
        self.Ra = - self.Rb - sum(P)
    
    
    def compute_shear_forces(self):
        T = []
        T.append(self.P[0])
        T.append(self.P[0] + self.Ra)
        T.append(T[1] + self.P[1])
        T.append(T[2] + self.P[2])
        T.append(T[3] + self.Rb)
        T.append(T[4] + self.P[3])
        self.T = T
        if T[4] + self.P[3] == 0.0:
            print('check ok')


    def compute_bending_moments(self):
        M = [0]
        M.append(self.a * self.T[0])
        M.append(M[1] + self.b[0] * self.T[1])
        M.append(M[2] + self.b[1] * self.T[2])
        M.append(M[3] + self.b[2] * self.T[3])
        M.append(round((M[4] + self.c * self.T[4]), 3))
        self.M = M
