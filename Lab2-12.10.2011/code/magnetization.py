import numpy as np

from matplotlib import pyplot as plt
import mc

#from scipy import constants

## temperature list
## from 1.9 to 4.0 with step 0.01
T_list = np.arange(1.5, 3.5, 0.01)

lattice_size = 25
sq_size = float(lattice_size ** 2)
eq_time = 100 ## time (in sweeps) when systems goes into equilibrium
meassure_time = eq_time/20 ## avoid autocorrelation
avg_number = None

sim = mc.MC(size = lattice_size, sweeps=5000, T = T_list[0], eq_time = eq_time)

magnetization_avg = []
energy_avg = []

for T in T_list:
    print T
    sim.T = T
    results = sim.run()

    #magnetization = map(lambda x: x / float(lattice_size**2), results['magnetization'])
    magnetization = results['magnetization'][0::meassure_time]
    if not avg_number:
        avg_number = float(len(magnetization) * sq_size)

    #file = open("data/M-T_%.2f.csv" % T, "w+")
    #file.writelines(["%d\n" % m for m in magnetization ])
    #file.close()

    #energy = results['energy'][eq_time:]

    ## compute monte carlo average magnetization for given temperature
    magnetization_avg.append(abs(sum(magnetization))/avg_number)

    #plt.plot([ m/avg_number for m in magnetization ])
    #plt.show()
    #energy_avg.append(sum(energy)/float(len(energy)*sq_size))


plt.title('Magnetization for N = %d' % lattice_size)
plt.ylabel('<M>')
plt.xlabel('T')
plt.plot(T_list, magnetization_avg)

plt.show()
