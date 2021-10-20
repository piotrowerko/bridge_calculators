#dobór zbrojenia teowego przekroju żelbetowego
#poszukiwanie granicznej wysokości strefy ściskanej
Hdz = 1.5
Bdz = 1
Lteo = 30
Hpl = 0.25
a = 1
MoblmaxToT = 1200
lambda_ = 0.8
fyd = 435 * 1000
eta = 1.0 #współczynnik równy jeden dla wybranej stali
fcd = 35.0 * 1000 / 1.4
epsilon_cu = 3.5 * 10 ** -3
epsilon_yd=2.1*10**-3  #fyd/Ey (~435/205)=2.1
ksi_eff_lim=lambda_*epsilon_cu/(epsilon_cu+epsilon_yd)
print('ksi_eff_lim=',ksi_eff_lim)
d = Hdz-0.035-0.012-0.032-0.035/2
print('d =',d)
x_eff_lim=ksi_eff_lim*d
print('x_eff_lim =', x_eff_lim)
As1 = eta * fcd * Bdz * x_eff_lim / fyd
print('As1 =',As1)
#liczba fi 32
#liczba_fi_32=As1/(3.14*(0.032**2)/4)
#print('liczba_fi_32=',liczba_fi_32)
#eta*fcd*Bdz*x_eff_lim*(d-0.5*x_eff_lim)
fcd

#szerekość efektywna przekroju:
b_eff=min(Lteo*0.2+Bdz,Hpl*12+Bdz,a)
print(Lteo*0.2+Bdz,Hpl*12+Bdz,a,'b_eff=',b_eff)

#ocena czy przekrój jest teowy czy pozornie teowy:
#obliczenie nośności przekroju w przypadku ściskania całej półki:
M_rd_p_eff=eta*fcd*b_eff*Hpl*(d-0.5*Hpl)
print('M_rd_p_eff=',M_rd_p_eff)
Beta=Hpl/d
M_rd_p_eff_=eta*fcd*b_eff*d*d*Beta*(1-0.5*Beta)
print('M_rd_p_eff_=',M_rd_p_eff_)
if MoblmaxToT<M_rd_p_eff:
    print('Przekrój pozornie teowy !')
else:
    print('Przekrój rzeczywiście teowy !')
    
mi_eff=MoblmaxToT/(fcd*b_eff*d**2)
print(mi_eff)
ksi_eff=0.03 #z tablicy 4.8
if ksi_eff<ksi_eff_lim:
    print('przekrój nie wymaga zbrojenia w strefie ściskanej')
else:
    print('trzeba liczyć też As2 !')

As1=eta*ksi_eff*d*b_eff*fcd/fyd
print('As1=',As1)

#liczba fi 32
liczba_fi_32=As1/(3.14*(0.032**2)/4)
print('liczba_fi_32=',liczba_fi_32)

procent_zbr=100*As1/(Bdz*Hdz)
print('procent_zbr=',procent_zbr,'%')