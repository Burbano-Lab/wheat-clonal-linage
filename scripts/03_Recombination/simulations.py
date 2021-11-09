import matplotlib.pyplot as plt
import numpy as np
import FFPopSim as h
import itertools

def simul(params):
    L = params['L']
    N = params['N']
    ### Set up
    pop = h.haploid_highd(L)
    pop.set_wildtype(N)
    pop.mutation_rate = params['mu']
    pop.outcrossing_rate = params['r']
    pop.crossover_rate = params['cr']
    # Half of the population with genotype: 00000... and half with genotype: 111111...
    pop.set_genotypes([[0]*L, [1]*L],[N*0.5, N*0.5])
    # Evolve population for a certain number of generations
    pop.evolve(params['gen'])
    # Compute pairwise LD
    LDs = []
    for pair in itertools.combinations(range(L), 2):
        l1 = pair[0]
        l2 = pair[1]
        diff = l2 - l1
        LDs.append([diff, abs(pop.get_LD(l1, l2))])
    LDs = np.array(LDs)
    # Compute median value per bin
    out = []
    for i in range(1,int(max(LDs[:,0]))+1):
        med = np.median(LDs[LDs[:,0] == i, 1])
        out.append([med])
    out = (LDs, np.array(out))
    return(out)

# PARAMETERS
#   gen: Number of generations
#   N: Population size
#   L: Number of loci
#   mu: mutation rate per site per generation
#   r: Probability of sexual reproduction per generation
#   cr: Crossover rate per site per generation
param = [{'gen':100, 'N':100000, 'L':300, 'mu':0.0003 , 'r':0.0, 'cr':0.01},
        {'gen':100, 'N':100000, 'L':300, 'mu':0.0003 , 'r':0.001, 'cr':0.01},
        {'gen':100, 'N':100000, 'L':300, 'mu':0.0003 , 'r':0.01, 'cr':0.01},
        {'gen':100, 'N':100000, 'L':300, 'mu':0.0003 , 'r':0.1, 'cr':0.01},
        {'gen':100, 'N':100000, 'L':300, 'mu':0.0003 , 'r':1.0, 'cr':0.01}]

# Run function and produce plots
cols = ['#ab62c0', '#6ea659', '#cb566b', '#648ace', '#c1823b', 'b']
show = 'both' # 'raw' / 'average' / 'both'
indx = np.array([[i] for i in range(1,param[0]['L'])])
fig = plt.figure()
ax = fig.add_subplot(1,1,1)
for p in range(len(param)):
    result_i = simul(param[p])
    if show == 'average':
        ax.plot(indx, result_i[1], color = cols[p], ls='-', lw=3, label=r''+str(param[p])+'')
    else:
        ax.scatter(result_i[0][:,0], result_i[0][:,1], color = cols[p], alpha=0.1, marker='.')
        if show == 'both':
            ax.plot(indx, result_i[1], color = cols[p], ls='-', lw=3, label=r''+str(param[p])+'')
plt.legend()
ax.set_ylim([0,0.25])
plt.xlabel('Distance between loci')
plt.ylabel('LD')
plt.ioff()
plt.show()
