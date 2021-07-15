import matplotlib.pyplot as plt
from belka import Belka

def main():
    belka_1 = Belka(2, [3, 5, 3], 4)
    belka_1.compute_reactions([-2, -5, -8, -3])
    belka_1.compute_shear_forces()
    belka_1.compute_bending_moments()
    print(belka_1.b)
    print(belka_1.P)
    print(f'reakcja Ra: {belka_1.Ra}')
    print(f'reakcja Rb: {belka_1.Rb}')
    print(f'siły tnące: {belka_1.T[:-1]}')
    print(f'momenty: {belka_1.M}')
    fig, ax = plt.subplots()
    ax.plot([0, belka_1.a, belka_1.a+ belka_1.b[0], belka_1.a+ belka_1.b[0]+belka_1.b[1],
             belka_1.a + sum(belka_1.b), belka_1.a + sum(belka_1.b) + belka_1.c], 
            [-belka_1.M[i] for i in range(len(belka_1.M))], label='bending moment')
    ax.axhline(linewidth=2, color='#d62728')
    ax.set_xlabel('Beam Length [m]')
    ax.set_ylabel('Bending moment [kNm]')
    plt.grid()
    plt.legend()
    plt.show()
    

if __name__ == '__main__':
    main()