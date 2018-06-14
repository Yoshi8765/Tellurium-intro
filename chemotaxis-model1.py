import tellurium as te
import numpy as np
import matplotlib.pyplot as plt
import roadrunner
roadrunner.Config.setValue(roadrunner.Config.MAX_OUTPUT_ROWS,1000000)

r = te.loada('''

# ligand binding and unbinding
J0: X -> C ;   a1*S*X
J1: C -> X ;   d1*C
J2: Xp -> Cp ; a1*S*Xp
J3: Cp -> Xp ; d1*Cp

# receptor methylation
J4: X -> Xp ;  a2*R*X/(X+KM2)
J5: Xp -> X ;  d2*B*Xp/(Xp+KD2)
J6: C -> Cp ;  a2*R*C/(C+KM2)
J7: Cp -> C ;  d2low*B*Cp/(Cp+KD2)

# rate constants
a1 = 500             # ligand association rate
d1 = 1               # ligand dissociation rate

# so when we have a sudden addition of chemoattractant

a2 = 0.01            # maximal methylation velocity
KM2 = 10^(-3)        # MM constant for methylation
d2 = 0.1             # maximal demethylation velocity of free receptor
d2low = 0.02         # maximal demethylation velocity of complexed receptor (< d2)
KD2 = 1              # MM constant for demethylation

vx = 1               # receptor activity in free state
vc = 0.1             # receptor activity in bound state

# initial conditions
X = 4                # total concentration of receptor complex (MCP/CheW/CheA)
C = 0
Xp = 0
Cp = 0
S = 0

R = 1               # total concentration of methylase (CheR)
B = 1            # total concentration of demethylase (CheB)

# step wise increase 
at time >= 400: S = 100
''')

#r.simulate(0,0.1,1000)
#r.plot()

tmax = 1500;

m = r.simulate(0,tmax,tmax)
T = m[:,0]
Xp = m[:,3]
Cp = m[:,4]

activity = r.vx*Xp + r.vc*Cp
plt.plot(T,activity)
plt.xlim(0,tmax)
plt.xlabel('time')
plt.ylabel('Total activity')

plt.show()