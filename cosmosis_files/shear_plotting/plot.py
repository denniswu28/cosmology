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
dir_name = "CosmoDC2_GW_01_12_test"
full_dir_name = "CosmoDC2"
r_dir = "/hpc/group/cosmology/denniswu/cosmosis_files/output/"+full_dir_name+"/data_vector_full/2pt_angle.txt"
r_cut_dir = "/hpc/group/cosmology/denniswu/cosmosis_files/output/"+dir_name+"/data_vector/2pt_angle.txt"
xi_data_dir = "/hpc/group/cosmology/denniswu/cosmosis_files/output/"+full_dir_name+"/data_vector_full/2pt_data.txt"
xi_data_cut_dir = "/hpc/group/cosmology/denniswu/cosmosis_files/output/"+dir_name+"/data_vector/2pt_data.txt"
xi_theory_dir = "/hpc/group/cosmology/denniswu/cosmosis_files/output/"+full_dir_name+"/data_vector_full/2pt_theory.txt"
xi_theory_cut_dir = "/hpc/group/cosmology/denniswu/cosmosis_files/output/"+dir_name+"/data_vector/2pt_theory.txt"
xi_bin1_dir = "/hpc/group/cosmology/denniswu/cosmosis_files/output/"+full_dir_name+"/data_vector_full/2pt_bin1.txt"
xi_bin2_dir = "/hpc/group/cosmology/denniswu/cosmosis_files/output/"+full_dir_name+"/data_vector_full/2pt_bin2.txt"
xi_bin1_cut_dir = "/hpc/group/cosmology/denniswu/cosmosis_files/output/"+dir_name+"/data_vector/2pt_bin1.txt"
xi_bin2_cut_dir = "/hpc/group/cosmology/denniswu/cosmosis_files/output/"+dir_name+"/data_vector/2pt_bin2.txt"
# sig_dir = "/hpc/group/cosmology/denniswu/cosmosis_files/output/CosmoDC2/data_vector/2pt_covariance.txt"

r = 60*180/3.1415*np.loadtxt(r_dir, comments="#", delimiter=",", unpack=False).flatten()
r_cut = 60*180/3.1415*np.loadtxt(r_cut_dir, comments="#", delimiter=",", unpack=False).flatten()
xi_data = np.loadtxt(xi_data_dir, comments="#", delimiter=",", unpack=False).flatten()
xi_data_cut = np.loadtxt(xi_data_cut_dir, comments="#", delimiter=",", unpack=False).flatten()
xi_theory = np.loadtxt(xi_theory_dir, comments="#", delimiter=",", unpack=False).flatten()
xi_theory_cut = np.loadtxt(xi_theory_cut_dir, comments="#", delimiter=",", unpack=False).flatten()
xi_bin1 = np.loadtxt(xi_bin1_dir, comments="#", delimiter=",", unpack=False).flatten()
xi_bin2 = np.loadtxt(xi_bin2_dir, comments="#", delimiter=",", unpack=False).flatten()
xi_bin1_cut = np.loadtxt(xi_bin1_cut_dir, comments="#", delimiter=",", unpack=False).flatten()
xi_bin2_cut = np.loadtxt(xi_bin2_cut_dir, comments="#", delimiter=",", unpack=False).flatten()

data_t = np.transpose(np.vstack((r,xi_data,xi_theory,xi_bin1,xi_bin2)))
data = np.vstack((r,xi_data,xi_theory,xi_bin1,xi_bin2))

data_t_cut = np.transpose(np.vstack((r_cut,xi_data_cut,xi_theory_cut,xi_bin1_cut,xi_bin2_cut)))
data_cut = np.vstack((r_cut,xi_data_cut,xi_theory_cut,xi_bin1_cut,xi_bin2_cut))

data_1_1 = np.transpose(data_t[np.logical_and(data[3] == 1, data[4] == 1)])
data_1_1_cut = np.transpose(data_t_cut[np.logical_and(data_cut[3] == 1, data_cut[4] == 1)])
# data_1_2 = np.transpose(data_t[np.logical_and(data[5] == 1, data[6] == 2)])
# data_1_3 = np.transpose(data_t[np.logical_and(data[5] == 1, data[6] == 3)])
# data_1_4 = np.transpose(data_t[np.logical_and(data[5] == 1, data[6] == 4)])
# data_2_2 = np.transpose(data_t[np.logical_and(data[5] == 2, data[6] == 2)])
# data_2_3 = np.transpose(data_t[np.logical_and(data[5] == 2, data[6] == 3)])
# data_2_4 = np.transpose(data_t[np.logical_and(data[5] == 2, data[6] == 4)])
# data_3_3 = np.transpose(data_t[np.logical_and(data[5] == 3, data[6] == 3)])
# data_3_4 = np.transpose(data_t[np.logical_and(data[5] == 3, data[6] == 4)])
# data_4_4 = np.transpose(data_t[np.logical_and(data[5] == 4, data[6] == 4)])
print(data_1_1)
print(data_1_1_cut)

# 13 5 7 15
# 10 4 7 15



fig, axs = plt.subplots(2, 2, figsize=(12, 12))
axs[0, 0].plot(r[0:20], xi_data[0:20], linestyle="",marker=".", color='black', label = r'CosmoDC2 Truth $\xi_+ (\theta)$')
# axs[0, 0].plot(r_cut[0:10], xi_data_cut[0:10], linestyle="",marker=".", color='red', label = r'CosmoDC2 after-cut Truth $\xi_+ (\theta)$')
axs[0, 0].plot(r[0:20], xi_theory[0:20], linestyle="solid",marker="", color='navy', label = r'CosmoDC2 best-fit theory $\xi_+ (\theta)$')
# axs[0, 0].plot(r_cut[0:10], xi_theory_cut[0:10], linestyle="solid",marker="", color='red', label = r'CosmoDC2 after-cut best-fit theory $\xi_+ (\theta)$')
axs[0, 0].set_yscale('log')
axs[0, 0].set_ylabel(r'$\xi_+$')
axs[0, 0].legend(loc = 'lower left')
axs[0, 1].plot(r[20:40], xi_data[20:40], linestyle="",marker=".", color='black', label = r'CosmoDC2 Truth $\xi_- (\theta)$')
# axs[0, 1].plot(r_cut[10:18], xi_data_cut[10:18], linestyle="",marker=".", color='red', label = r'CosmoDC2 after-cut Truth $\xi_- (\theta)$')
axs[0, 1].plot(r[20:40], xi_theory[20:40], linestyle="solid",marker="", color='navy', label = r'CosmoDC2 best-fit theory $\xi_- (\theta)$')
# axs[0, 1].plot(r_cut[10:18], xi_theory_cut[10:18], linestyle="solid",marker="", color='red', label = r'CosmoDC2 after-cut best-fit theory $\xi_- (\theta)$')
axs[0, 1].set_yscale('log')
axs[0, 0].set_ylabel(r'$\xi_-$')
axs[0, 1].legend(loc = 'lower left')
axs[1, 0].plot(r[40:60], xi_data[40:60], linestyle="",marker=".", color='black', label = r'CosmoDC2 Truth $\gamma_t (\theta)$')
axs[1, 0].plot(r_cut[0:7], xi_data_cut[0:7], linestyle="",marker=".", color='red', label = r'CosmoDC2 after-cut Truth $\gamma_t (\theta)$')
axs[1, 0].plot(r[40:60], xi_theory[40:60], linestyle="solid",marker="", color='navy', label = r'CosmoDC2 best-fit theory $\gamma_t (\theta)$')
axs[1, 0].plot(r_cut[0:7], xi_theory_cut[0:7], linestyle="solid",marker="", color='red', label = r'CosmoDC2 after_cut best-fit theory $\gamma_t (\theta)$')
axs[1, 0].set_yscale('log')
axs[0, 0].set_ylabel(r'$\gamma_t$')
axs[1, 0].legend(loc = 'lower left')
axs[1, 1].plot(r[60:80], xi_data[60:80], linestyle="",marker=".", color='black', label = r'CosmoDC2 Truth $\omega (\theta)$')
axs[1, 1].plot(r_cut[7:19], xi_data_cut[7:19], linestyle="",marker=".", color='red', label = r'CosmoDC2 after-cut Truth $\omega (\theta)$')
axs[1, 1].plot(r[60:80], xi_theory[60:80], linestyle="solid",marker="", color='navy', label = r'CosmoDC2 best-fit theory $\omega (\theta)$')
axs[1, 1].plot(r_cut[7:19], xi_theory_cut[7:19], linestyle="solid",marker="", color='red', label = r'CosmoDC2 after-cut best-fit theory $\omega (\theta)$')
axs[1, 1].set_yscale('log')
axs[0, 0].set_ylabel(r'$w_\theta$')
axs[1, 1].legend(loc = 'lower left')
for ax in axs.flat:
    ax.set(xlabel='Angle (arcmin)')
    ax.set_xscale('log')
fig.suptitle('3x2pt correlations')
fig.savefig('/hpc/group/cosmology/denniswu/cosmosis_files/shear_plotting/shear_2.png')

