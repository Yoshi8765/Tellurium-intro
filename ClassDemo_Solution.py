import numpy as np
import matplotlib.pyplot as plt
import tellurium as te
import roadrunner
import antimony

r = te.loada("""
    J0: -> M ; a_m
    J1: M -> ; d_m*M
    J2: -> P ; a_p*M
    J3: P -> ; d_p*P
    
    # Parameters
    a_m = 10; d_m = 1
    a_p = 500; d_p = 0.05
    
    # Initial values
    M = 0
    P = 0
""")

result = r.simulate(0,200,1000) #(start, end, timepoints)

# Find the steady-state value 
P_ss = result[-1,2]
# Find the response time (time to 1/2 steady state)
count = 0 # result index tracker
for i in result[:,2]:
    if i >= 0.5*P_ss:
        t_r = result[count,0]
        break
    else:
        count += 1

r.plot() # quick plot

# Plotting using pyplot
plt.figure(1)
plt.subplot(211)
plt.plot(result[:,0], result[:,1], 'r', label = "mRNA")
plt.ylabel('Copies', fontsize='12') 
plt.xlabel('Time (hr)', fontsize='12')
plt.legend(loc=4)

plt.subplot(212)
plt.plot(result[:,0], result[:,2], 'b', label = "Protein")
plt.plot([0,200], [P_ss,P_ss], 'k--', label = "Steady-State")
plt.plot([0,t_r], [0.5*P_ss,0.5*P_ss], 'k')
plt.plot([t_r,t_r], [0, 0.5*P_ss], 'k', label = "Response Time")
plt.ylim([0,120000])
plt.ylabel('Copies', fontsize='12') 
plt.xlabel('Time (hr)', fontsize='12')
plt.legend(loc=4)
plt.tight_layout()

#%% Doing a parametric sweep using a for loop

# Investigate the input output relationship between 
# steady state Protein levels (P_ss)
# and mRNA production rate (a_m)
# or mRNA decay rate (d_m)

a_m = np.linspace(1,50,50)
d_m = np.linspace(0.02, 1,50)
P_ss = np.zeros(len(d_m))
T_r = np.zeros(len(d_m))
count = 0
for i in d_m:
    r = te.loada("""
        J0: -> M ; a_m
        J1: M -> ; d_m*M
        J2: -> P ; a_p*M
        J3: P -> ; d_p*P
        
        # Parameters
        a_m = 10;
        a_p = 500; d_p = 0.05
        
        # Initial values
        M = 0
        P = 0
        """
        # Iterating parameters
        "d_m = " + str(i) + ";"
    )
    result = r.simulate(0,200,1000)
    
    # Find the steady-state value
    P_ss[count] = result[-1,2]
    
    # Find the response time
    index = 0
    for j in result[:,2]:
        if j >= 0.5*result[-1,2]:
            t_r = result[index,0]
        break
    else:
        index += 1
    T_r[count] = t_r # Store the response time
    count += 1
    print("d_m value:" + str(count))

#%% Plots
P_ss2 = r.a_m*r.a_p/r.d_p/d_m

plt.figure(2)
plt.plot(d_m,P_ss,label="Numerical")
plt.plot(d_m,P_ss2,label="Analytical")
plt.ylabel('Steady State Value', fontsize='12') 
plt.xlabel('mRNA decay rate', fontsize='12')
plt.legend()
    
plt.figure(3)
plt.plot(d_m,T_r,label="Response Time")
plt.ylabel('Response Time', fontsize='12') 
plt.xlabel('mRNA decay rate', fontsize='12')
