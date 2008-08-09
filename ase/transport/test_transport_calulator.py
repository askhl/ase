from ase.transport.calculators import TransportCalculator
from ase.transport.tools import tri2full
import numpy as npy


H_lead = npy.zeros([4,4])

# On-site energies are zero
for i in range(4):
    H_lead[i,i] = 0.0

# Nearest neighbor hopping is -1.0
for i in range(3):
    H_lead[i,i+1] = -1.0
    H_lead[i+1,i] = -1.0

# Next-nearest neighbor hopping is 0.2
for i in range(2):
    H_lead[i,i+2] = 0.2
    H_lead[i+2,i] = 0.2

H_scat = npy.zeros([6,6])
# Principal layers on either side of S
H_scat[:2,:2] = H_lead[:2,:2]
H_scat[-2:,-2:] = H_lead[:2,:2]

# Scattering region
H_scat[2,2] = 0.0
H_scat[3,3] = 0.0
H_scat[2,3] = -0.8
H_scat[3,2] = -0.8

# External coupling
H_scat[1,2] = 0.2
H_scat[2,1] = 0.2
H_scat[3,4] = 0.2
H_scat[4,3] = 0.2

pl = 2
energies = npy.arange(-3,3,0.02)
#print "creating"
#calc = TransportCalculator()
#print "setting pl"
#calc.set(pl=pl)
#print "setting h" 
#calc.set(h=H_scat)
#print "setting h1"
#calc.set(h1=H_lead)
#print "setting h2"
#calc.set(h2=H_lead)
#print "setting energies"
#calc.set(energies=energies)
calc = TransportCalculator(pl=pl,
                           h=H_scat,
                           h1=H_lead,
                           h2=H_lead,
                           energies=energies)
#print "setting pdos"
calc.trans.set(pdos=[0,1])
#save the original hamiltonian before working on it
h = calc.h_pp.copy()
s = calc.s_pp.copy()
#"print calculating T"
T = calc.get_transmission()
dos = calc.get_dos()
pdos = calc.get_pdos()

#subdiagonalize
ha, sa, u, eps = calc.subdiagonalize_bfs([0,1])
calc.set(h=ha,s=sa)
Ta = calc.get_transmission()
dosa = calc.get_dos()
pdosa = calc.get_pdos()

#remove coupling
hb, sb = calc.cutcoupling_bfs([0])
calc.set(h=hb,s=sb)
Tb = calc.get_transmission()
dosb = calc.get_dos()
pdosb = calc.get_pdos()

#revert to the original h,s
calc.set(h=h,s=s)
Tc = calc.get_transmission()
dosc = calc.get_dos()
pdosc = calc.get_pdos()

#import pylab
#pylab.plot(energies,T)
#pylab.plot()
