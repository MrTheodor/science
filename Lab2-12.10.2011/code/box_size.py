
from matplotlib import pyplot as plt
import mc

#from scipy import constants

# box size
Box_size = range(5, 21)

T_c = 2.269

eq_time = 1000 ## time (in sweeps) when systems goes into equilibrium
meassure_time = eq_time/20 ## avoid autocorrelation

magnetization_avg = []
energy_avg = []

for lattice_size in Box_size:
    print lattice_size
    sq_size = float(lattice_size ** 2)
    sim = mc.MC(size = lattice_size, sweeps=4000, T = T_c, eq_time = eq_time)

    results = sim.run()

    #magnetization = map(lambda x: x / float(lattice_size**2), results['magnetization'])
    magnetization = results['magnetization'][0::meassure_time]
    
    
    #file = open("data/M-T_%.2f.csv" % T, "w+")
    #file.writelines(["%d\n" % m for m in magnetization ])
    #file.close()

    #energy = results['energy'][eq_time:]

    ## compute monte carlo average magnetization for given temperature
    magnetization_avg.append(abs(sum(magnetization))/float(len(magnetization)*sq_size))
    #energy_avg.append(sum(energy)/float(len(energy)*sq_size))
    del sim


plt.title('Magnetization for T = %.2f' % T_c)
plt.ylabel('<M>')
plt.xlabel('Lattice size')
plt.plot(Box_size, magnetization_avg)

plt.show()
