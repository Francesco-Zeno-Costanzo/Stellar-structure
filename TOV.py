import numpy as np
import matplotlib.pyplot as plt
from scipy.interpolate import interp1d

#costanti in CGS
G = 6.6730831e-8
c = 2.99792458e10
Msun = 1.988435e33

#fattori di conversione vari ed eventuali
Ba_to_Pa = 1/10
cm_to_km = 1/100/1000

MeV_fm3_to_ba = 1.6021766e33
MeV_fm3_to_g_cm3 = MeV_fm3_to_ba/c**2


class TOV:
    """
    Classe per la risoluzione della tov;
    le quantita in input DEVONO avere dimenioni di:
    densità di energia[MeV/fm^3], pressione[MeV/fm^3].
    Il codice converte le quantità nel sistema CGS per fare la simulazione.
    La massa e il raggio della stella vengono restituiti
    in unità di masse solari e km rispettivamente.
    Per pressione e densità è possibile scegliere se CGS o SI.
    """

    def __init__(self, ene, pre):
        """
        Converte le quantità da [MeV/fm^3]
        a [g/cm^3] per la densità di energia
        a [g/(cms^2)] per la pressione
        Vine poi eseguita un'interpolazione cubica

        Parametri
        ---------
        ene : 1d array
            densità di energia in [MeV/fm^3]
        pre : 1d array
            pressione in [MeV/fm^3]
        """
        # converto in unità cgs
        ene *= MeV_fm3_to_g_cm3
        pre *= MeV_fm3_to_ba

        #interpolo per avere la densità energia in funzione della pressione
        self.rho = interp1d(pre, ene, kind='cubic')
        #interpolo per avere la pressione in funzione della densità di energia
        self.press = interp1d(ene, pre, kind='cubic')

        #eos letta da file in unità cgs
        self.__rho = ene
        self.__pre = pre

        #range in cui data la eos si possono fare simulazioni
        self.min_rho = np.min(ene)
        self.max_rho = np.max(ene)

        self.min_p = np.min(pre)
        self.max_p = np.max(pre)


    def plot_eos(self):
        """
        Esegue il plot dell'equazione di stato letta
        da file e dell'interpolazione eseguita
        """
        conv1 = MeV_fm3_to_g_cm3
        conv2 = MeV_fm3_to_ba

        #Divido per i fattori di conversione in modo da
        #riavere unità di misura presenti nel file inziale

        fig1 = plt.figure(1, figsize=(15,7))
        frame1=fig1.add_axes((.1,.35,.4,.6))
        #plot eos
        plt.plot(self.__rho/conv1, self.__pre/conv2, marker='o', label='EOS')
        #plot dell'interoplazione
        plt.plot(self.__rho/conv1, self.press(self.__rho)/conv2, marker='^', label='EOS_int')
        plt.ylabel(r'$P [MeV/fm^{3}]$')
        plt.legend()
        plt.grid()

        frame2=fig1.add_axes((.55,.35,.4,.6))
        #plot eos
        plt.plot(self.__pre/conv2, self.__rho/conv1, marker='o', label='EOS')
        #plot dell'interoplazione
        plt.plot(self.__pre/conv2, self.rho(self.__pre)/conv1, marker='^', label='EOS_int')
        plt.ylabel(r'$\varepsilon [MeV/fm^{3}]$')
        plt.legend()
        plt.grid()

        #per controllare abbia interpolato bene
        frame3=fig1.add_axes((.1,.1,.4,.2))
        #differanza fra eos ed eos interpolata
        plt.plot(self.__rho/conv1, (self.press(self.__rho)-self.__pre)/conv2, marker='.', linestyle='')
        plt.xlabel(r'$\varepsilon [MeV/fm^{3}]$')
        plt.ylabel('EOS-EOS_int')
        plt.grid()

        frame4=fig1.add_axes((.55,.1,.4,.2))
        #differanza fra eos ed eos interpolata
        plt.plot(self.__pre/conv2, (self.rho(self.__pre)-self.__rho)/conv1, marker='.', linestyle='')
        plt.xlabel(r'$P [MeV/fm^{3}]$')
        plt.ylabel('EOS-EOS_int')
        plt.grid()

        plt.show()


    def tov_eq(self, r, P, m):
        """
        Equazione da risolvere

        Parametri
        ---------
        r : float
            poszione all'interno della stella
        P : float
            Pressione della stella a raggio r
        M : float
            Massa della stella a raggio r

        Returns
        ---------
        dPdr : float
            derivata della pressione calcolata in r
        dmdr : float
            derivata della massa calcolata in r
        """
        #posso simulare solo nel range di validità della eos
        if P < self.min_p or P > self.max_p:
            return 0., 0.

        rho = self.rho(P)

        dPdr = -G*(rho + P/c**2.)*(m + 4.0*np.pi*r**3*P/c**2.)
        dPdr = dPdr/(r*(r - 2.0*G*m/c**2.))

        dmdr = 4.0*np.pi*rho*r**2.

        return dPdr, dmdr


    def controllo(self, dens):
        """
        Controlla se la densità centrale iserita è nel
        range di validita dell'equazione di stato.
        Se non lo è viene sollevata l'eccezione
        restituendo il messaggio di errore.

        Parametri
        ---------
        dens : float
            densità in [g/cm^3]
        """
        #posso simulare solo nel range di validità della eos
        conv = MeV_fm3_to_g_cm3
        if dens < self.min_rho or dens > self.max_rho:
            a = f'La densità centrale scelta: {dens/conv:.5e} è fuori dal rand ella EOS.\n'
            b = f'Il vaore deve essere compreso tra:{self.min_rho/conv:.5e} e {self.max_rho/conv:.5e}'
            errore = f'{a}{b}'
            raise Exception(errore)


    def RK4(self, f, P_c, M_c, dr):
        """
        Calcola la soluzione usando runge_kutta4.

        Parametri
        ---------
        f : callable
            funzione da integrare.
        P_c : float
            pressione centrale.
        M_c : float
            Massa centrale.
        dr : float
            passo dell'integratore.

        Returns
        ---------
        r : 1d array
            coordinata radiale nella stella in cm.
        P : 1d array
            pressione nella stella in funzione di r in Ba.
        M : 1d array
            massa nella stella in funzione di r in g.
        """

        P = []
        M = []
        r = []

        P.append(P_c)
        M.append(M_c)
        r.append(dr)

        while P[-1] > self.min_p:
            pk1, mk1 = f(r[-1], P[-1], M[-1])
            pk2, mk2 = f(r[-1] + dr/2, P[-1] + pk1*dr/2, M[-1] + mk1*dr/2)
            pk3, mk3 = f(r[-1] + dr/2, P[-1] + pk2*dr/2, M[-1] + mk2*dr/2)
            pk4, mk4 = f(r[-1] + dr, P[-1] + pk3*dr, M[-1] + mk3*dr)

            P_n = P[-1] + (dr/6)*(pk1 + 2*pk2 + 2*pk3 + pk4)
            M_n = M[-1] + (dr/6)*(mk1 + 2*mk2 + 2*mk3 + mk4)
            r_n = r[-1] + dr

            P.append(P_n)
            M.append(M_n)
            r.append(r_n)

        r = np.array(r)
        P = np.array(P)
        M = np.array(M)

        return r, P, M


    def solve(self, dens_c, dr, SI=False):
        """
        Calcola le codizioni inziali su P e M e chiama RK4

        Parametri
        ---------
        dens_c : float
            densità centrale in [MeV/fm^3]
        dr : float
            passo di integrazione
        SI : bool, oprional
            se True la pressione viene restituita in unità SI

        Return
        ---------
        r : 1d array
            coordinata radiale nella stella in km.
        P : 1d array
            pressione nella stella in funzione di r in Ba/Pa .
        M : 1d array
            massa nella stella in funzione di r in Masse solari.
        """

        dens_c *= MeV_fm3_to_g_cm3
        #controllo se possibile eseguire la simulazione
        self.controllo(dens_c)
        #condizioni iniziali
        p = self.press(dens_c)
        rho = self.rho(p)
        m = (4/3)*np.pi*rho*dr**3.

        r, P, M = self.RK4(self.tov_eq, p, m, dr)

        r = r*cm_to_km
        M = M/Msun

        if SI:
            P = P*Ba_to_Pa

        return r, P, M