# -*- coding: utf-8 -*-
"""
Spyder Editor

This is a temporary script file.
"""

import tellurium as te
import numpy
import matplotlib.pyplot as plt
# close all previous plots
plt.close('all')
# Loading Tellurium, differential equations, and initializing variables.
r = te.loada ('''

    x1' = a*x1 + b*x2;
    x2' = c*x1 + d*x2;
    
   x1 = .5; x2 = .5; a = -4; b = -2; c = -2; d = 4;
    
''')
# User input for matrix coefficients & assigned to existing Tellurium variables
r.a = float(raw_input("what is the a coefficient? "))
r.b = float(raw_input("what is the b coefficient? "))
r.c = float(raw_input("what is the c coefficient? "))
r.d = float(raw_input("what is the d coefficient? "))

# creating mmesh grid and defining diferential equations for streamplot call
Y,X=numpy.mgrid[-6:6:100j, -6:6:100j]
u=r.a*X+r.b*Y
v=r.c*X+r.d*Y
# using roadrunner to simulate/solve the differential equations a set tikme interval 
# set number of datapoints.
m = (r.simulate (0, 15, 600, ['x1', 'x2'])).T
# creating time vector for time course plot
t=numpy.linspace(0,15,600)
# jacobian constructed
jacobian = numpy.array([[r.a,r.b],[r.c,r.d]])
# eigen values determined and rounded to 10 decimal places
eigen, vectors = numpy.linalg.eig(jacobian)
eigen[0]=numpy.around(eigen[0], decimals=10, out=None)
eigen[1]=numpy.around(eigen[1], decimals=10, out=None)
# Printing for user to see jacobian and eigen values
print "The jacobian is: "
print jacobian
print ""
print "The eigen values are: "
print eigen
print ""
# determination of behavior based on the real and imaginary parts of eigen values
# a string is initialized with the behavior which is used for stating behavior in
# console and plot titles
behavior = str()
if eigen[0].imag == 0 or eigen[1].imag == 0:
    if eigen[0].real > 0 and eigen[1].real > 0:
        behavior = "Unstable Behavior"        
        print "This system will have an " + behavior + "."
elif eigen[0].real > 0 or eigen[1].real > 0:
        behavior = "Saddle Behavior"
        print "This system will have a " + behavior + "."
elif eigen[0].real < 0 and eigen[1].real < 0:
        behavior = "Stable Behavior"
        print "This system will have a " + behavior + "."
           
elif eigen[0].real < 0 and eigen[1].real < 0:
    behavior = "Spiral Sink Behavior"    
    print "This system will have a " + behavior + "."
elif  eigen[0].real > 0 or eigen[1].real > 0:
    behavior = "Spiral Source Behavior"    
    print "This system will have a " + behavior + "."
else:
     behavior = "Center Behavior"
     print "This system will have a " + behavior + "."


# plotting Argan plot and making it look good
plt.scatter(eigen.real,eigen.imag)
plt.xlabel('Real Part')
plt.ylabel('Imaginary Part')
plt.legend()
plt.title('Argand Plot - ' + behavior)
plt.grid(True, which='both')
plt.show()
#creating a subplot for time course plot and linear phase portrait
f,ax = plt.subplots(1,2)
# plotting time course plot
ax[1].plot(t, m[0, :], label="x1")
ax[1].plot(t, m[1, :], label="x2")
ax[1].set_xlabel('Time (s)')
ax[1].set_ylabel('Concentration')
ax[1].legend()
ax[1].set_title('Time Course Plot - ' + behavior)
ax[1].grid(True, which='both')
# plotting linear phase plot
ax[0].streamplot(X,Y,u,v,density=[1,1])
ax[0].set_xlabel('x1')
ax[0].set_ylabel('x2')
ax[0].legend()
ax[0].set_title('Linear Phase Portrait - ' + behavior)
ax[0].axis([-6,6,-6,6])
ax[0].grid(True, which='both')
