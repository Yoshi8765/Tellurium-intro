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
J4: C -> Cp ;  a2*R*C/(C+KM2)

# receptor demethylation
J5: Xp -> X ;  d2*B*Xp/(Xp+KD2)

# rate constants
a1 = 100
d1 = 100

a2 = 0.004
KM2 = 10^(-8)
KD2 = 1
d2 = 0.6

# initial conditions
X = 4
C = 0
Xp = 0
Cp = 0
S = 1

B = 1       # total concentration of methylase (CheR)
R = 1      # total concentration of demethylase (CheB)

# events
at time >= 400: S = 100
''')

#r.simulate(0,0.1,1000)
#r.plot()

tmax = 1500;
vx = 10;


m = r.simulate(0,tmax,tmax)
T = m[:,0]
Xp = m[:,3]

activity = vx*Xp


plt.plot(T,activity)
plt.xlim(0,tmax)
plt.xlabel('time')
plt.ylabel('activity')


plt.show()