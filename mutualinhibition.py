import tellurium as te
import numpy as np
import matplotlib.pyplot as plt

# model conditions

r = te.loada('''
J0:   -> x; (v1*s)/((1+s)*(1 + 4*y + 6*w1*pow(y,2) + 4*pow(w1,2)*pow(y,3) + pow(w1,3)*pow(y,4)))
J1: x ->  ; d1*x
J2:   -> y; (v2*s)/((1+s)*(1 + 4*x + 6*w2*pow(x,2) + 4*pow(w2,2)*pow(x,3) + pow(w2,3)*pow(x,4)))
J3: y -> ;  d2*y

# kinetic constants
s = 1      # the concentration of signal present
v1 = 1         # maximal rate of x expression
v2 = 1         # maximal rate of y expression
w1 = 10        # cooperativity in inhibition of x by y
w2 = 10        # cooperativity in inhibition of y by x
d1 = 1         # degradation rate of x
d2 = 1         # degradation rate of y

#initial conditions
x = 0.05
y = 0.1
''')

# nullcline plots
R = r.simulate(0,100,1000)

t = R[:,0]
x = R[:,1]
y = R[:,2]

# time traces
plt.plot(t,x, color="green", label="x")
plt.plot(t,y, color="blue", label="y")
plt.legend(loc="upper left")
plt.hold('ON')
plt.ylim([0, 1])
plt.xlabel('time')
plt.ylabel('levels')
plt.show()

# compute expression for x-nullcline (dx/dt = 0)
yn1 = np.linspace(0,1,200)
xn1 = (r.v1*r.s)/((1+r.s)*(1 + 4.*yn1 + 6.*r.w1*pow(yn1,2) + 4.*pow(r.w1,2)*pow(yn1,3) + pow(r.w1,3)*pow(yn1,4)))/(r.d1)

# compute expression for y-nullcline (dy/dt = 0)
xn2 = np.linspace(0,1,200)
yn2 = (r.v2*r.s)/((1+r.s)*(1 + 4.*xn2 + 6.*r.w2*pow(xn2,2) + 4.*pow(r.w2,2)*pow(xn2,3) + pow(r.w2,3)*pow(xn2,4)))/(r.d2)


# phase plane analysis
plt.hold('ON')
plt.plot(x,y,color="black",label="trajectory")
plt.plot(xn1,yn1,color="green", label="x-nullcline")
plt.plot(xn2,yn2,color="blue", label="y-nullcline")
plt.legend(loc="upper left")

#plt.xlim([0, 0.5])
#plt.ylim([0, 0.5])
plt.xlabel('x (levels)')
plt.ylabel('y (levels)')
plt.show()
