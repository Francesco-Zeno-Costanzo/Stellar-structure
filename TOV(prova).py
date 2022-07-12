import numpy as np
import matplotlib.pyplot as plt
from scipy import constants as cst

## CAVEAT:
# Il codice risolve la tov allo stesso modo delle equazioni di lane - Emden
# si è assunto c = G = 1 e non si sa come come tornare ad unità fisiche,
# in più anche i valori dei parametri q o P_c sono scelti a caso non sapendo
# la fisica per determinarli nelle "unità di misura" adottate nel codice;
# gli andamenti delle soluzioni sembrano però corretti se confrontati con
# quanto si trova in rete da simulazioni di articoli seri.
##

#qui c'è ancora parvenza di fisica
hbar = cst.value(u'reduced Planck constant')
mn = cst.value(u'neutron mass')
K = ((3*np.pi)**(2/3)*hbar**2)/(5*mn**(8/3))
n = 1.5

#da qui in poi si entra nella tana del bianconiglio
q = 0.5


def TOV(r, p, m):
    '''
    equazioni da risolvere con dipendenza della
    politropica già  inserita
    '''
    dPdr = -((q + p)/(n+1))*(m + 4*np.pi*r**3*p**(n+1)/q)/(r*(r-2*m))

    dmdr = 4*np.pi*r**2*p**n

    return dPdr, dmdr

def RK4(P_c, M_c, x_c):
    '''
    soluzione tramite il metodo runge-kutta di ordine 4
    '''

    dx = 1e-4

    P = []
    M = []
    x = []

    P.append(P_c)
    M.append(M_c)
    x.append(x_c)

    #pressione "zero" identifica la fine della stella
    while P[-1] > 1e-5:
        pk1, mk1 = TOV(x[-1], P[-1], M[-1])
        pk2, mk2 = TOV(x[-1] + dx/2, P[-1] + pk1*dx/2, M[-1] + mk1*dx/2)
        pk3, mk3 = TOV(x[-1] + dx/2, P[-1] + pk2*dx/2, M[-1] + mk2*dx/2)
        pk4, mk4 = TOV(x[-1] + dx, P[-1] + pk3*dx, M[-1] + mk3*dx)

        P_n = P[-1] + (dx/6)*(pk1 + 2*pk2 + 2*pk3 + pk4)
        M_n = M[-1] + (dx/6)*(mk1 + 2*mk2 + 2*mk3 + mk4)
        x_n = x[-1] + dx

        #P potrebbe diventare negativa durante
        #il passo e quindi la derivta comlessa
        #causa sostituzione della politropica
        if type(P_n) == complex: break

        P.append(P_n)
        M.append(M_n)
        x.append(x_n)

    return x, P, M


#singola stella
x, p, m = RK4(1, 0, 1e-5)

P = []
M = []
X = []

for ra, ma , pa in zip(x, m, p):
    md = ma/np.sqrt((q*K)**(-n))
    rd = ra/np.sqrt((q*K)**(-n))
    pd = (pa**(n+1))*(K**(-n)*q**(-n-1))
    M.append(md)
    X.append(rd)
    P.append(pd)

plt.figure(1)
plt.suptitle('massa e pressione in funzione del raggio')
plt.subplot(121)
plt.grid()
plt.xlabel('$r$')
plt.ylabel('$M$')
plt.plot(X, M)
plt.subplot(122)
plt.grid()
plt.plot(X, P)
plt.xlabel('$r$')
plt.ylabel('pressione')
plt.show()

##Curva massa raggio

R = []
M = []
for i in np.arange(0.001, 0.02, 0.0005):
    x, p, m = RK4(i, 0, 1e-5)

    md = m[-1]/np.sqrt((q*K)**(-n))
    rd = x[-1]/np.sqrt((q*K)**(-n))

    M.append(md)
    R.append(rd)


for i in np.arange(0.02, 0.35, 0.01):
    x, p, m = RK4(i, 0, 1e-5)

    md = m[-1]/np.sqrt((q*K)**(-n))
    rd = x[-1]/np.sqrt((q*K)**(-n))

    M.append(md)
    R.append(rd)


plt.figure(2)
plt.title('Massa in funzione del raggio')
plt.xlabel('R')
plt.ylabel('M')
plt.grid()
plt.plot(R, M, marker='.')

plt.show()
