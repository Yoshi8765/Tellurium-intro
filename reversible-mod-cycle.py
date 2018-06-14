import tellurium as te
import numpy as np
import matplotlib.pyplot as plt

r = te.loada('''
J0: X -> Xp; VM0*S*X/(X+KM0)
J1: Xp-> X ; VM1*E1*Xp/(Xp+KM1)

# kinetic constants
S = 10  # this is E0
E1 = 10
VM0 = 1
VM1 = 1.5
KM0 = 3
KM1 = 5

# initial conditions
X  = 10 
Xp = 0
''')

# k is independent variable
k_values = []
# steady state values of A
Xp_values = []
for S in np.linspace(0,50,num=500):
    r.resetAll()
    r.S = S
    # call simulate first to approach correct fixed point
    r.simulate(0,100)
    # now call steady state
    #r.steadyState()
    k_values.append(S)
    Xp_values.append(r.Xp)
    
# plotting code
plt.plot(k_values, Xp_values, marker='.')
plt.ylim([-0.1, 10.1])
plt.xlim([0, 20])
plt.tick_params(labelsize=14)
plt.xlabel('S',fontsize=20)
plt.ylabel('Xp (steady state)',fontsize=20)
plt.show()