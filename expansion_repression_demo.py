import tellurium as te
import numpy as np
import matplotlib.pyplot as plt


# lattice points 
N = 50


# simulation timepoints
timepoints = 50000
tmax = 500000

# we will fill in these terms later
antimonyStr = '''
  # M diffusion terms
  {Mdiff}
  # M degradation terms
  {Mdeg}

  # initial values
  {Minitial}
  {Einitial}

  # E diffusion terms
  {Ediff}
  # E production terms
  {Eprod}
  # E degradation terms
  {Edeg}

 # global values
  L = 200 # length
  N = {N} # number of lattice points
  dx = L/N # length of lattice cell
  
  # morphogen parameters
  Dm = 10 # diffusion constant
  dm = Dm/dx^2 # diffusion rate constant
  Sm =  20  #source production rate
  am  = 0.1 # M degradation rate constant 
  Ke = 0.1    # expander concentration required for inhibition of morphogen degradation
  
  # expander parameters
  De = 10 # diffusion constant   
  de= De/dx^2 # diffusion rate constant
  Km = 0.02 # morphogen concentration required for repression of expander
  h = 4 # hill coefficient for expander repression
  be = 0.01    #   E production rate constant
  ae = 0.00001  # E degradation rate constant  
'''

# expressions to be substituted
blocks = {'N': N}

######################### MORPHOGEN TERMS ###################
# build the diffusion terms
diff_terms = []
# add first term
diff_terms.append('M_diff0: -> c_0; dm*(Sm + c_1 - c_0)')
for i in range(1,N-1):
    diff_terms.append('M_diff{i}: -> c_{i}; dm*(c_{inext} + c_{iprev} - 2*c_{i})'.format(i=i, inext=i+1, iprev=i-1))
# add the last term
diff_terms.append('M_diff{last}: -> c_{last}; dm*(c_{prev} - c_{last})'.format(last=N-1, prev=N-2))
blocks['Mdiff'] = '\n  '.join(diff_terms)

# build the morphogen degradation terms
deg_terms = []
for i in range(N):
    deg_terms.append('M_deg{i}: c_{i} ->; (Ke*am*c_{i})/(Ke+e_{i})'.format(i=i))       # morphogen 
blocks['Mdeg'] = '\n  '.join(deg_terms)

# build the initial values
initial_terms = []
for i in range(N):
    initial_terms.append('c_{i} = 0'.format(i=i))
blocks['Minitial'] = '\n  '.join(initial_terms)

######################## EXPANDER TERMS #####################
# build the diffusion terms
diff_terms = []
# add first term
diff_terms.append('E_diff0: -> e_0; de*(e_1 - e_0)')
for i in range(1,N-1):
    diff_terms.append('E_diff{i}: -> e_{i}; de*(e_{inext} + e_{iprev} - 2*e_{i})'.format(i=i, inext=i+1, iprev=i-1))
# add the last term
diff_terms.append('E_diff{last}: -> e_{last}; dm*(e_{prev} - e_{last})'.format(last=N-1, prev=N-2))
blocks['Ediff'] = '\n  '.join(diff_terms)

# build the expander production terms
prod_terms = []
for i in range(N):
    prod_terms.append('E_prod{i}: -> e_{i}; be*(Km^h/(Km^h+c_{i}^h))'.format(i=i))       # morphogen repression
blocks['Eprod'] = '\n  '.join(prod_terms)

# build the expander degradation terms
deg_terms = []
for i in range(N):
    deg_terms.append('E_deg{i}: e_{i} ->; ae*e_{i}'.format(i=i))
blocks['Edeg'] = '\n  '.join(deg_terms)

# build the initial values
initial_terms = []
for i in range(N):
    initial_terms.append('e_{i} = 0.1'.format(i=i))
blocks['Einitial'] = '\n  '.join(initial_terms)


######################### SIMULATION ##########################

print(antimonyStr.format(**blocks))
m = te.loada(antimonyStr.format(**blocks))

r = m.simulate(0,tmax,timepoints)

# plot the results
ind = np.linspace(1,timepoints-1,3)
x = np.linspace(0,m.L,m.N)

# plot the 
for i in ind:   # choose all the time points for display
    j = int(i)   
    t1 = r[j,0]   # timepoint
    m = r[j,1:(N+1)]   #
    e = r[j,(N+1):2*N+1] 
    plt.plot(x,m)   # plot morphogen
    plt.ylabel('morphogen levels')
    plt.ylim([0.01, 300])
    plt.hold(True)
    plt.plot((0, max(x)),(1, 1))
    plt.yscale('log')  
    plt.title('time ='+str(t1))
    plt.show()
      
    plt.plot(x,e)   # plot the expander
    plt.xlabel('distance (x)')
    plt.ylabel('expander levels')
    plt.ylim([0, 1])
    plt.show()

    




