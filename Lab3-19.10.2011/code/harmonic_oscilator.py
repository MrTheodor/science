## Techniki symulacyjne
## Jakub Krajniak <jkrajniak@gmail.com>
##
## 19.10.2011
##

import math

from matplotlib import pyplot as plt
import numpy as np

dt = (2*math.pi) / 2000
omega_dt = dt
omega_dt_sq = omega_dt**2

x0 = 1.0
v0 = 0.001
x1 = x0 - v0*dt

K = 5000

sim_vec = [x0, x1]

for k in range(1, K ):
    value = float(2 - omega_dt_sq) * sim_vec[k] - sim_vec[k-1]
    sim_vec.append(value)

plt.plot(np.sin(np.arange(0, K+1)*(2*math.pi/2000)+0.5*math.pi))
plt.plot(sim_vec)

plt.legend(["Analytic", "MD"])

plt.title("MD vs analytic solution for model of harmonic oscilator")

plt.show()
