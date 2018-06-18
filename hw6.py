# -*- coding: utf-8 -*-
"""
Created on Tue Feb 21 22:23:27 2017

@author: Administrator
"""
import tellurium as te
import roadrunner
import numpy as np
from sys import exit
from scipy import signal
from matplotlib import pyplot as plt

# user inputs to make the jacobian
a=input('Enter A: ') 
b=input('Enter B: ') 
c=input('Enter C: ') 
d=input('Enter D: ') 
jac=np.matrix([[a,b], [c,d]])

# Set up the U and V domain meshgrid
X,Y=np.mgrid[-4:4:400j,-4:4:400j]
X1=jac[0,0]*X+jac[0,1]*Y
X2=jac[1,0]*X+jac[1,1]*Y

eigw,eigv=np.linalg.eig(jac)
#using tellurium to model time dependent plots
r = te.loada(
'x\'='+str(jac[0,0])+'*x+'+str(jac[0,1])+'*y;\n'+ #dx1/dt
'y\'='+str(jac[1,0])+'*x+'+str(jac[1,1])+'*y;\n'+ #dx2/dt
'x=4;y=-4;') #Initial Conditions
m=r.simulate(0,8,400) #running the model

#making a phase plot with streamlines
plt.subplots(1,2, figsize=(10,4))
plt.subplot(121)
plt.xlabel('x',fontsize='16')
plt.ylabel('y',fontsize='16')
plt.streamplot(Y,X,X2,X1,density=[1,1])
plt.ylim((-4,4))
plt.xlim((-4,4))
#determine system stability type
if eigw[0].imag == 0 and eigw[1].imag == 0:
    if eigw[0] > 0 and eigw[1] > 0:
        print "System Stability: Unstable System."
    elif eigw[0] < 0 and eigw[1] < 0:
        print "System Stability: Stable system."
    elif eigw[0] == 0 and eigw[1] == 0:
        print "System Stability: Constant valued system."
    else:
        print "System Stability: Saddle Point exists."
elif eigw[0].real == 0 and eigw[1].real == 0:
    print "System Stability: Center exists."
else:
    if eigw[0].real > 0 and eigw[1].real > 0:
        print "System Stability: Unstable spiral system."
    if eigw[0].real < 0 and eigw[1].real < 0:
        print "System Stability: Stable spiral system."

exit('lol')
# Scan over zeta, a parameter for a second-order system
zetaRange = np.arange(0.1,1.1,0.1)
#{
f1 = plt.figure()
for i in range(0,9):
    den = [1, 2*zetaRange[i], 1]
    print den
    s1 = signal.lti(num, den)
    # Specify our own frequency range: np.arange(0.1, 5, 0.01)
    w, mag, phase = signal.bode(s1, np.arange(0.1, 5, 0.01).tolist())
    plt.semilogx (w, mag, color="blue", linewidth="1")
plt.xlabel ("Frequency")
plt.ylabel ("Magnitude")
plt.savefig ("c:\\mag.png", dpi=300, format="png")

plt.figure()

for i in range(0,9):
    den = [1, zetaRange[i], 1]
    s1 = signal.lti(num, den)
    w, mag, phase = signal.bode(s1, np.arange(0.1, 10, 0.02).tolist())
    plt.semilogx (w, phase, color="red", linewidth="1.1")
plt.xlabel ("Frequency")
plt.ylabel ("Amplitude")
plt.savefig ("c:\\phase.png", dpi=300, format="png")


r = te.loada("""
x1'=jac[1]*x1+jac[2]*x2;
x2'=jac[3]*x1+jac[4]*x2;
x1=[]
""")
