"""
This script is to compute the 2-point shear-shear correlation of the threshholded truth catalog using treecorr. The lines reading in the catalog is adapted from Prof. Troxel's script.

Author: Dennis Wu
"""

# All the libs needed
# import fitsio as fio
# import csv
# import twopoint
# import numpy as np
# import matplotlib
# matplotlib.use ('agg')
# import matplotlib.pyplot as plt
# import galsim
# import pylab
# import healpy as hp
# import os
# import treecorr
# import glob
# import re
# import utilities as util
# import math


import pandas as pd
import numpy as np
from scipy import stats
from sklearn.model_selection import LeaveOneOut
from sklearn.metrics import accuracy_score
import matplotlib
import matplotlib.pyplot as plt
matplotlib.use ('agg')


def load_data(stars_path, gals_path):
    stars = pd.read_csv(stars_path)
    gals = pd.read_csv(gals_path)
    return stars, gals


stars_path = '/hpc/group/cosmology/denniswu/scolnic/stars.csv'
gals_path = '/hpc/group/cosmology/denniswu/scolnic/gals.csv'

stars, gals = load_data(stars_path, gals_path)


print("max z", np.max(stars['z']))

plt.hist(stars['z'], label = r'stars', density = True, bins = 50)
plt.hist(gals['z'], label = r'gals', density = True, bins = 50)
plt.xlabel(r'$z$')
plt.ylabel(r'$dN/dz$')
plt.legend(loc = 'upper right')
plt.title("Redshift dist of galaxies samples and QSO samples")
plt.savefig('/hpc/group/cosmology/denniswu/scolnic/gal_star_hist.pdf')
plt.clf()

# obj_type = 'roman'
# twopoint_file_path = '/hpc/group/cosmology/denniswu/cosmosis_files/data/pseudo_'+obj_type+'.fits'
# plot_out_dir =  '/hpc/group/cosmology/denniswu/twopoint_plots/'+obj_type
# T1 = twopoint.TwoPointFile.from_fits(twopoint_file_path)
# T1.plots(plot_out_dir, latex = False)



# fr = ['Y106','J129','H158','F184']
# dec_min = -42
# dec_max = -38
# ra_min = 51
# ra_max = 56

# obj_type = 'truth'
# lens_type = 'len_'+obj_type
# source_type = 'source_'+obj_type
# bin_num = 5

# gg_file_path = '/hpc/group/cosmology/denniswu/data/gg_corr_all.fits'
# ng_file_path = '/hpc/group/cosmology/denniswu/data/ng_corr_21.fits'
# ng_rand_file_path = '/hpc/group/cosmology/denniswu/data/ng_rand_21.fits'
# nn_file_path = '/hpc/group/cosmology/denniswu/data/nn_corr_21.fits'
# cov_file_path = '/hpc/group/cosmology/denniswu/covmat/truth/cov_truth'
# full_file_path = '/hpc/group/cosmology/denniswu/data/pseudo.fits'
# full_cov_file_path = '/hpc/group/cosmology/denniswu/cosmosis_files/data/pseudo.fits'

# bins_tomo = [1,2,3,4,5]
# nbins_tomo = len(bins_tomo)

# nbins_nz = 149

# nbins = 20
# min_sep, max_sep = 1, 100
# bin_size = math.log(max_sep / min_sep)/ nbins        
# logr = np.linspace(0, nbins*bin_size, nbins, endpoint=False, dtype=float)
# logr += math.log(min_sep) + 0.5*bin_size
    
# rnom = np.exp(logr)
# half_bin = np.exp(0.5*bin_size)
# left_edges = rnom / half_bin
# right_edges = rnom * half_bin


# ''' Write lens and source nofz to twopoint fits file. '''
# # Fetch nz dist from truth, create hist
# nz_dir = '/hpc/group/cosmology/denniswu/hist_file.npy'
# nz_source = np.load(nz_dir,allow_pickle='TRUE').item()[source_type]
# nz_lens = np.load(nz_dir,allow_pickle='TRUE').item()[lens_type]

# hist_source, bins_source = zip(*nz_source)
# hist_lens, bins_lens = zip(*nz_lens)
# bins_source = bins_source[0]
# bins_lens = bins_lens[0]

# n_source = [np.sum(arr) for arr in list(hist_source)]
# n_lens = [np.sum(arr) for arr in list(hist_lens)]

# bins_source_mid = [(bins_source[i]+bins_source[i+1])/2 for i in range(nbins_nz)]
# bins_lens_mid = [(bins_lens[i]+bins_lens[i+1])/2 for i in range(nbins_nz)]

# print(len(np.repeat(bins_tomo, [x * nbins for x in reversed(range(1, nbins_tomo + 1))])))
# print(np.repeat(bins_tomo, [x * nbins for x in reversed(range(1, nbins_tomo + 1))]))

# print(len(np.repeat(np.concatenate([bins_tomo[i:] for i in range(nbins_tomo)]), nbins)))
# print(np.repeat(np.concatenate([bins_tomo[i:] for i in range(nbins_tomo)]), nbins))


# obj_type = 'roman'
# lens_type = 'len_'+obj_type
# source_type = 'source_'+obj_type

# nz_dir = '/hpc/group/cosmology/denniswu/hist_file.npy'

# # file = np.load(nz_dir,allow_pickle='TRUE').item()
# # print(file.keys()) ##to list the available keys 

# nz_source = np.load(nz_dir,allow_pickle='TRUE').item()[source_type]
# nz_lens = np.load(nz_dir,allow_pickle='TRUE').item()[lens_type]

# hist_source, bins_source = zip(*nz_source)
# hist_lens, bins_lens = zip(*nz_lens)

# bins_source = bins_source[0]
# bins_lens = bins_lens[0]

# nbins_nz = len(bins_source)-1

# bins_source_mid = [(bins_source[i]+bins_source[i+1])/2 for i in range(nbins_nz)]
# bins_lens_mid = [(bins_lens[i]+bins_lens[i+1])/2 for i in range(nbins_nz)]

# n_source = np.divide([np.sum(arr) for arr in list(hist_source)],(20*3600))
# n_lens = np.divide([np.sum(arr) for arr in list(hist_lens)],(20*3600))

# print(n_source)
# print(n_lens)

# source_file = np.column_stack((np.array(bins_source[:-1]), *hist_source))
# lens_file = np.column_stack((np.array(bins_source[:-1]), *hist_source))

# # Save to a file
# np.savetxt("nz/source_"+ obj_type + ".nz", source_file, fmt='%.18e', delimiter=' ')
# np.savetxt("nz/lens_"+ obj_type + ".nz", lens_file, fmt='%.18e', delimiter=' ')






# def extract_coordinates_regex(directory_string):
#     # Regular expression to find the pattern before '.fits.gz'
#     match = re.search(r'(\d+\.\d+_-?\d+\.\d+)\.fits\.gz$', directory_string)
#     if match:
#         return match.group(1)
#     return None

# def extract_unique_entries(data):
#     dtype_obj = data[0].dtype
#     unique_entries = {}
#     for obj in data:
#         if obj['x']>0 and obj['gal_star'] == 0 and np.logical_and(np.logical_and(obj['ra'] < ra_max, obj['ra'] > ra_min),np.logical_and(obj['dec'] < dec_max, obj['dec'] > dec_min)):
#             coordinates = (obj['ra'], obj['dec'])  # Extract coordinates, assuming dictionary format
#             if coordinates not in unique_entries:
#                 unique_entries[coordinates] = obj  # Store the whole entry
#         else:
#             continue
#     return np.array(list(unique_entries.values()), dtype = dtype_obj)

# dc2_coaddlist = fio.FITS('/hpc/group/cosmology/phy-lsst/public/dc2_sim_output/truth/coadd/dc2_coaddlist.fits.gz')[-1].read()
# # print(dc2_coaddlist.dtype)
# # print(dc2_coaddlist[1].dtype)
# # for row in dc2_coaddlist:
# #     print(row.shape)
# #     print(row['coadd_i'])
# #     break
# dc2_coaddlist_dict = {row[0]: row for row in dc2_coaddlist}


# for i,f in enumerate(np.sort(glob.glob('/hpc/group/cosmology/phy-lsst/public/dc2_sim_output/truth/coadd/dc2_index_*.fits.gz'))):
#     try:
#         print(f)
#         tmp = fio.FITS(f)[-1].read(columns=['ra','dec', 'mag_J129','mag_F184','mag_H158','mag_Y106', 'ind', 'gal_star','x'])
#         coordinate = extract_coordinates_regex(f)
#         unique = extract_unique_entries(tmp)
#         print(coordinate)
#         mask = dc2_coaddlist_dict[coordinate]
#         # print(mask)
#         # print(len(tmp))
#         # print(max(tmp['ra']), mask['coadd_ra']+mask['d_ra'], min(tmp['ra']), mask['coadd_ra']-mask['d_ra'])
#         # print(max(tmp['dec']), mask['coadd_dec']+mask['d_dec'], min(tmp['dec']), mask['coadd_dec']-mask['d_dec'])
#         tmp = tmp[np.logical_and(np.logical_and(tmp['ra'] < mask['coadd_ra']+mask['d_ra']/2, tmp['ra'] > mask['coadd_ra']-mask['d_ra']/2),
#                                     np.logical_and(tmp['dec'] < mask['coadd_dec']+mask['d_dec']/2, tmp['dec'] > mask['coadd_dec']-mask['d_dec']/2))]
#         print(len(tmp))
#         print(len(unique))
        
#     except:
#         pass
# print(tmp[1])
# print(unique[1])

# dc2_truth_names = np.sort(glob.glob('/hpc/group/cosmology/phy-lsst/public/dc2_sim_output/truth/coadd/dc2_index_*.fits.gz'))
# for i in range(10):
#     print(dc2_truth_names[i])
#     coordinate = extract_coordinates_regex()
#     print(dc2_coaddlist_dict[coordinate])
#     mask = dc2_coaddlist_dict[extract_coordinates_regex(f)]


# # Read in and cut truth objects
# dc2_truth = util.read_truth_catalog('/hpc/group/cosmology/phy-lsst/public/dc2_sim_output/truth/coadd/dc2_index_*.fits.gz')
# dc2_truth = dc2_truth[dc2_truth['x']>0] # select valid
# dc2_truth = dc2_truth[dc2_truth['gal_star'] == 0]
# dc2_truth = dc2_truth[np.logical_and(np.logical_and(dc2_truth['ra'] < ra_max, dc2_truth['ra'] > ra_min),np.logical_and(dc2_truth['dec'] < dec_max, dc2_truth['dec'] > dec_min))]

# dc2_truth_old= util.read_truth_catalog_old('/hpc/group/cosmology/phy-lsst/public/dc2_sim_output/truth/coadd/dc2_index_*.fits.gz')
# dc2_truth_old = dc2_truth_old[dc2_truth_old['x']>0] # select valid
# dc2_truth_old = dc2_truth_old[dc2_truth_old['gal_star'] == 0]
# dc2_truth_old = dc2_truth_old[np.logical_and(np.logical_and(dc2_truth_old['ra'] < ra_max, dc2_truth_old['ra'] > ra_min),np.logical_and(dc2_truth_old['dec'] < dec_max, dc2_truth_old['dec'] > dec_min))]

# print(dc2_truth.size)
# print(dc2_truth_old.size)

# # Read in and cut dectection objects
# dc2_det = util.read_det_catalog('/hpc/group/cosmology/phy-lsst/public/dc2_sim_output/detection/dc2_det_*.fits.gz')
# dc2_det = dc2_det[dc2_det['number']>0] # select valid
# dc2_det = dc2_det[(dc2_det['flags']<4)&(dc2_det['flux_auto']/dc2_det['fluxerr_auto']>5)] # select flag=0 and flux/err ratio >5
# dc2_det = dc2_det[(dc2_det['alphawin_j2000']>ra_min)&(dc2_det['alphawin_j2000']<ra_max)&(dc2_det['deltawin_j2000']>dec_min)&(dc2_det['deltawin_j2000']<dec_max)]

# # Match both catalogs
# dc2_det_match, dc2_truth_match, mask = util.get_match(dc2_det, dc2_truth)
# print(dc2_truth_match.dtype.names)
# print(dc2_det_match.dtype.names)
# print(len(dc2_truth_match))
# print(len(dc2_det_match))
# print(dc2_truth_match['ra'][0], dc2_truth_match['dec'][0])
# print(dc2_det_match['alphawin_j2000'][0], dc2_det_match['deltawin_j2000'][0])
# print(dc2_truth_match['ra'][1], dc2_truth_match['dec'][1])
# print(dc2_det_match['alphawin_j2000'][1], dc2_det_match['deltawin_j2000'][1])