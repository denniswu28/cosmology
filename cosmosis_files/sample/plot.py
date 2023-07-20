"""
This script is to compute the 2-point count-count correlation of the threshholded truth catalog using treecorr with Landy-Szalay estimator. The lines reading in the catalog is adapted from Prof. Troxel's script.

Author: Dennis Wu
"""

# All the libs needed
import numpy as np
import matplotlib
matplotlib.use ('agg')
import matplotlib.pyplot as plt

# plot the corr functions for NN
r_dir = "/hpc/group/cosmology/denniswu/phy555/output/des-default/data_vector/2pt_angle.txt"
xi_data_dir = "/hpc/group/cosmology/denniswu/phy555/output/des-default/data_vector/2pt_data.txt"
xi_theory_dir = "/hpc/group/cosmology/denniswu/phy555/output/des-default/data_vector/2pt_theory.txt"
xi_theory_1_dir = "/hpc/group/cosmology/denniswu/phy555/output/des-1/data_vector/2pt_theory.txt"
xi_theory_2_dir = "/hpc/group/cosmology/denniswu/phy555/output/des-2/data_vector/2pt_theory.txt"
xi_bin1_dir = "/hpc/group/cosmology/denniswu/phy555/output/des-default/data_vector/2pt_bin1.txt"
xi_bin2_dir = "/hpc/group/cosmology/denniswu/phy555/output/des-default/data_vector/2pt_bin2.txt"
sig_dir = "/hpc/group/cosmology/denniswu/phy555/output/des-default/data_vector/pantheon_covariance.txt"

r = 60*180/3.1415*np.loadtxt(r_dir, comments="#", delimiter=",", unpack=False).flatten()
xi_data = np.loadtxt(xi_data_dir, comments="#", delimiter=",", unpack=False).flatten()
xi_theory = np.loadtxt(xi_theory_dir, comments="#", delimiter=",", unpack=False).flatten()
xi_theory_1 = np.loadtxt(xi_theory_1_dir, comments="#", delimiter=",", unpack=False).flatten()
xi_theory_2 = np.loadtxt(xi_theory_2_dir, comments="#", delimiter=",", unpack=False).flatten()
xi_bin1 = np.loadtxt(xi_bin1_dir, comments="#", delimiter=",", unpack=False).flatten()
xi_bin2 = np.loadtxt(xi_bin2_dir, comments="#", delimiter=",", unpack=False).flatten()

data_t = np.transpose(np.vstack((r,xi_data,xi_theory,xi_theory_1,xi_theory_2,xi_bin1,xi_bin2)))
data = np.vstack((r,xi_data,xi_theory,xi_theory_1,xi_theory_2,xi_bin1,xi_bin2))

# binned_data = np.zero(4,4)

data_1_1 = np.transpose(data_t[np.logical_and(data[5] == 1, data[6] == 1)])
data_1_2 = np.transpose(data_t[np.logical_and(data[5] == 1, data[6] == 2)])
data_1_3 = np.transpose(data_t[np.logical_and(data[5] == 1, data[6] == 3)])
data_1_4 = np.transpose(data_t[np.logical_and(data[5] == 1, data[6] == 4)])
data_2_2 = np.transpose(data_t[np.logical_and(data[5] == 2, data[6] == 2)])
data_2_3 = np.transpose(data_t[np.logical_and(data[5] == 2, data[6] == 3)])
data_2_4 = np.transpose(data_t[np.logical_and(data[5] == 2, data[6] == 4)])
data_3_3 = np.transpose(data_t[np.logical_and(data[5] == 3, data[6] == 3)])
data_3_4 = np.transpose(data_t[np.logical_and(data[5] == 3, data[6] == 4)])
data_4_4 = np.transpose(data_t[np.logical_and(data[5] == 4, data[6] == 4)])
print(data_1_1)
n=np.linspace(0, 226, 227)

plt.plot(n, xi_data, linestyle="",marker=".", color='black', label = r'DES Real Data $\xi (\theta)$')
plt.plot(n, xi_theory, linestyle="solid",marker="", color='navy', label = r'DES LCDM theory $\xi (\theta)$')
plt.plot(n, xi_theory_1, linestyle="solid",marker="", color='red', label = r'DES $\Omega_m = 1$ theory $\xi (\theta)$')
plt.plot(n, xi_theory_2, linestyle="solid",marker="", color='green', label = r'DES $\Omega_k = 0.1$ theory $\xi (\theta)$')
# plt.plot(data_1_1[0], data_1_1[3], linestyle="solid",marker="", color='red', label = r'DES $\Omega_m = 1$ prediction $\xi (\theta)$')
# plt.plot(data_1_1[0], data_1_1[4], linestyle="solid",marker="", color='blue', label = r'DES $\Omega_k = 0.1$ prediction $\xi (\theta)$')



# plt.xscale('log')
plt.yscale('log')
plt.ylim(1e-10,1e-2)
# plt.xlabel(r'$\theta$ (arcmin)')
plt.xlabel('num')
plt.ylabel(r'$\xi$')
plt.legend(loc = 'lower left')
plt.title("2pt cosmic shear correlations")

plt.savefig('sample/shear.png')
plt.clf()
