import numpy as np
import matplotlib.pyplot as plt
from TOV import TOV


##EOS

path = r"C:\Users..."

#Unità di misura in MeV/fm^3, MeV/fm^3 e 1/fm^3
ene_dens, pre, num_dens = np.loadtxt(path, unpack=True)

tov = TOV(ene_dens, pre)
tov.plot_eos()

##singola stella

dr = 100
dens_c = 876
r, P, M = tov.solve(dens_c, dr, SI=True)


plt.figure(1, figsize=(10,6))
plt.suptitle(f'Simulazioe con densità centrale {dens_c} MeV/fm^{3}', fontsize=15)

plt.subplot(121)
plt.title('Massa in funzione del raggio')
plt.plot(r, M)
plt.ylabel('$M/M_{\odot}$')
plt.xlabel('r [km]')
plt.grid()

plt.subplot(122)
plt.title('Pressione in funzione del raggio')
plt.plot(r, P)
plt.xlabel('r [km]')
plt.ylabel('P [Pa]')
plt.grid()

plt.show()



##curva massa raggio

M = []
R = []

dr = 100

for dens_c in np.logspace(2.05, 3.2, 100):
#for dens_c in np.linspace(110, 1600, 80):
    r, _, m = tov.solve(dens_c, dr)
    M.append(m[-1])
    R.append(r[-1])


plt.figure(2, figsize=(9, 6))
plt.title('Curva Massa vs Raggio', fontsize=15)
plt.plot(R, M)
plt.grid()
plt.ylabel(r'$M/M_{\odot}$')
plt.xlabel(r'R [km]')


#salvo su file
path = r"C:\Users..."
file = open(path, "w")
riga = f'# Massa [M/Msun] \t Raggio [km]\n'
file.write(riga)
for m, r in zip(M, R):
    riga = f'{m} \t {r} \n'
    file.write(riga)
file.close()


plt.show()
