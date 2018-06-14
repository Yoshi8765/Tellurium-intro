import tellurium as te
import numpy as np
import matplotlib.pyplot as plt

w = 200   # cooperativity value
v1 = 21   # rate of transcription in the bound states
v0 = 1.5   # rate of transcription for the unbound states

amax = 0.12   # maximum on the x-axis of the plot

# Note here we changed variables a = K[A]
# to non-dimensionalize TF concentrations

# ploting the functional dependence of transcription rate on a
a = np.linspace(0,amax,5000)
v = v0 + (v1-v0)*(3*a + (2*w+1)*(a**2)+(w**2)*(a**3))/(1 + 3*a + (2*w+1)*(a**2) + (w**2)*(a**3))
plt.plot(a,v,marker='.')
plt.show()