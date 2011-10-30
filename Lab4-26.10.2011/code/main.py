from matplotlib import pyplot as plt
import mc
from scipy import constants

Temp = 2.0

sim = mc.MC(size=10, sweeps=50000, T = Temp, eq_time = 100)

results = sim.run()

magnetization = map(lambda x: x / float(sim.lattice_size_sq), results['magnetization'][0::150])
energy = map(lambda x: x / float(sim.lattice_size_sq), results['energy'][0::150])


plt.figure(1)
plt.subplot(211)
plt.title('Magnetization per spin for T=%f' % Temp)
plt.ylabel('<M>')
plt.plot(magnetization, '.')

plt.subplot(212)
plt.title('Energy per spin')
plt.xlabel('MC Sweeps')
plt.ylabel('Energy/NN')
plt.plot(energy)

plt.show()
