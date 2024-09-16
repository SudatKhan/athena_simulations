import numpy as np
from scipy.ndimage import shift
import sys
np.set_printoptions(threshold=sys.maxsize)
import matplotlib.pyplot as plt
sys.path.insert(0, '/home/sudat/athena/vis/python')
import athena_read

data = athena_read.athdf(f'/home/sudat/cca/single/disk.out1.00001.athdf')
r = data.get('x1v')

file = input(f"Enter directory path: ")
title = input('Enter Plot Title: ')
inner_sma = float(input(f'Enter the semimajor axis of the inner planet: '))
outer_sma = float(input(f'Enter the semimajor axis of the outer planet: '))
data1 = athena_read.athdf(file + f'/disk.out1.00001.athdf')
r1 = data1.get('x1v')

shifted_dens_sum = 0
for i in range(600,800,5):
    data = athena_read.athdf(f'/home/sudat/cca/single/disk.out1.00{i}.athdf')
    densp1 = data.get('dens')[0]
    phip1_c = data.get('x2f')
    time2p1 = data.get('Time')

    period1 = 2*np.pi*np.sqrt(1**3 /(1))
    planet_phi = 2*(np.pi / period1)*(time2p1)
    
    coeff = planet_phi/np.pi
    k = np.floor(coeff/2)
    beta = ((coeff * np.pi) - 2*np.pi*k)
    index = np.where(phip1_c > beta)[0][0]

    dens = shift(densp1, (-1*index, 0), mode="wrap")

    shifted_dens_sum += dens
shift_dens_avg= shifted_dens_sum / len(range(600,800,5))

shifted_dens_sum2 = 0
for i in range(600,800,5):
    data = athena_read.athdf(file + f'/disk.out1.00{i}.athdf')
    densp1 = data.get('dens')[0]
    phip1_c = data.get('x2f')
    time2p1 = data.get('Time')

    period1 = 2*np.pi*np.sqrt(inner_sma**3 /(1))
    planet_phi = 2*(np.pi / period1)*(time2p1)
    
    coeff = planet_phi/np.pi
    k = np.floor(coeff/2)
    beta = ((coeff * np.pi) - 2*np.pi*k)
    index = np.where(phip1_c > beta)[0][0]

    dens = shift(densp1, (-1*index, 0), mode="wrap")

    shifted_dens_sum2 += dens
shift_dens_avg2 = shifted_dens_sum2 / len(range(600,800,5))

shifted_dens_sum3 = 0
for i in range(600,800,5):
    data = athena_read.athdf(file + f'/disk.out1.00{i}.athdf')
    densp1 = data.get('dens')[0]
    phip1_c = data.get('x2f')
    time2p1 = data.get('Time')

    period1 = 2*np.pi*np.sqrt(outer_sma**3 /(1))
    planet_phi = 2*(np.pi / period1)*(time2p1)
    
    coeff = planet_phi/np.pi
    k = np.floor(coeff/2)
    beta = ((coeff * np.pi) - 2*np.pi*k)
    index = np.where(phip1_c > beta)[0][0]

    dens = shift(densp1, (-1*index, 0), mode="wrap")

    shifted_dens_sum3 += dens
shift_dens_avg3 = shifted_dens_sum3 / len(range(600,800,5))

plt.plot(r, np.std(shift_dens_avg, axis=0), label ='single')
plt.plot(r1, np.std(shift_dens_avg2, axis=0), label = 'planet 1')
plt.plot(r1, np.std(shift_dens_avg3, axis=0), label = 'planet 2')
plt.xscale('log')
plt.ylabel(r'$\sigma_{\Sigma}$')
plt.xlabel(f'Radius')
plt.title(title)
plt.legend()
plt.show()