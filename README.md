# Stellar-structure

Here there are codes for the resolution of the stellar structure in hydrostatic equilibrium in the framework the Newtonian gravity of general relativity.


## Newtonian gravity ##

Let's start with the Newtonian case, the equations are:

<img src="https://latex.codecogs.com/svg.image?\frac{dm(r)}{dr}&space;=&space;4&space;\pi&space;r^2&space;\rho(r)\hspace{5&space;mm}\frac{dP(r)}{dr}&space;=&space;\frac{-Gm(r)&space;\rho(r)}{r^2}" title="https://latex.codecogs.com/svg.image?\frac{dm(r)}{dr} = 4 \pi r^2 \rho(r)\hspace{5 mm}\frac{dP(r)}{dr} = \frac{-Gm(r) \rho(r)}{r^2}" />

Assuming a polytropic equation of state it is possible to obtain the Lane Emden equation:

<img src="https://latex.codecogs.com/svg.image?&space;\frac{1}{\xi^2}&space;\frac{d}{d\xi}&space;\left({\xi^2&space;\frac{d\theta}{d\xi}}\right)&space;=&space;-&space;\theta^n&space;" title="https://latex.codecogs.com/svg.image? \frac{1}{\xi^2} \frac{d}{d\xi} \left({\xi^2 \frac{d\theta}{d\xi}}\right) = - \theta^n " />

from which it is then possible to derive the limiting mass of white dwarfs, as is done in the code 'Chandrasekhar_limit.py'.


## General Relativity ##

In the context of general relativity the equations become instead:

<img src="https://latex.codecogs.com/svg.image?\\{}\frac{dm(r)}{dr}&space;=&space;4&space;\pi&space;r^2&space;\rho(r)\\\frac{dP(r)}{dr}&space;=&space;\frac{-Gm(r)&space;\rho(r)}{r^2}\Biggl(1&plus;&space;\frac{P(r)}{c^2&space;\rho(r)}\Biggr)\Biggl(1&plus;&space;\frac{4&space;\pi&space;r^3P(r)}{m(r)}\Biggr)\Biggl(1-&space;\frac{2m(r)G}{rc^2}\Biggr)^{-1}&space;" title="https://latex.codecogs.com/svg.image?\\{}\frac{dm(r)}{dr} = 4 \pi r^2 \rho(r)\\\frac{dP(r)}{dr} = \frac{-Gm(r) \rho(r)}{r^2}\Biggl(1+ \frac{P(r)}{c^2 \rho(r)}\Biggr)\Biggl(1+ \frac{4 \pi r^3P(r)}{m(r)}\Biggr)\Biggl(1- \frac{2m(r)G}{rc^2}\Biggr)^{-1} " />

In the code 'TOV(prova).py' these equations are solved as in the previous case for Lane-Emden with analogous calculations; and as before it is possible to calculate the limiting mass for, in this case, neutron stars.
Does the code work? Yes, the trends of the curves are correct and probably, perhaps, also the numerical values ​​(in some system of measurement); but having assumed an equation of state we lose generality of it. Let's try to fix it.
