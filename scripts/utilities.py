import fitsio as fio
import csv
import re
import numpy as np
import matplotlib
matplotlib.use ('agg')
import matplotlib.pyplot as plt
import treecorr
import glob
import scipy as sp
from scipy import spatial

def extract_coordinates_regex(directory_string):
    # Regular expression to find the pattern before '.fits.gz'
    match = re.search(r'(\d+\.\d+_-?\d+\.\d+)\.fits\.gz$', directory_string)
    if match:
        return match.group(1)
    return None

def extract_unique_entries(data):
    dtype_obj = data[0].dtype
    unique_entries = {}
    for obj in data:
        if obj['x']>0 and obj['gal_star'] == 0: #and np.logical_and(np.logical_and(obj['ra'] < ra_max, obj['ra'] > ra_min),np.logical_and(obj['dec'] < dec_max, obj['dec'] > dec_min)):
            coordinates = (obj['ra'], obj['dec'])  # Extract coordinates, assuming dictionary format
            if coordinates not in unique_entries:
                unique_entries[coordinates] = obj  # Store the whole entry
        else:
            continue
    return np.array(list(unique_entries.values()), dtype = dtype_obj)

def get_match(m,t):
    m_ori = m
    t_ori = t
    tree = spatial.cKDTree(np.c_[t['ra'].ravel(), t['dec'].ravel()])
    dd,ii = tree.query(np.c_[m['alphawin_j2000'].ravel(), m['deltawin_j2000'].ravel()], k=3)
    m1=dd[:,0]*60.*60.<1
    m2=dd[:,1]*60.*60.<1
    m3=dd[:,2]*60.*60.<1
    dm1=(t_ori['mag_F184'][ii[:,0]]-m_ori['mag_auto_F184'])
    dm2=(t_ori['mag_F184'][ii[:,1]]-m_ori['mag_auto_F184'])
    dm3=(t_ori['mag_F184'][ii[:,2]]-m_ori['mag_auto_F184'])
    mask = np.ones(len(m_ori))*-1
    mask[m1] = ii[m1,0]
    mask[m2&(np.abs(dm2)<np.abs(dm1))]= ii[m2&(np.abs(dm2)<np.abs(dm1)),1]
    mask[m3&(np.abs(dm3)<np.abs(dm1))&(np.abs(dm3)<np.abs(dm2))]= ii[m3&(np.abs(dm3)<np.abs(dm1))&(np.abs(dm3)<np.abs(dm2)),2]
    return m[mask>=0],t[np.unique(mask[mask>=0].astype(int))]

def get_match_old(m,t):
    m_ori = m
    t_ori = t
    tree = spatial.cKDTree(np.c_[m['alphawin_j2000'].ravel(), m['deltawin_j2000'].ravel()])
    dd,ii = tree.query(np.c_[t['ra'].ravel(), t['dec'].ravel()], k=3)
    m1=dd[:,0]*60.*60.<1
    m2=dd[:,1]*60.*60.<1
    m3=dd[:,2]*60.*60.<1
    dm1=(m_ori['mag_auto_F184'][ii[:,0]]-t_ori['mag_F184'])
    dm2=(m_ori['mag_auto_F184'][ii[:,1]]-t_ori['mag_F184'])
    dm3=(m_ori['mag_auto_F184'][ii[:,2]]-t_ori['mag_F184'])
    mask = np.ones(len(t_ori))*-1
    mask[m1] = ii[m1,0]
    mask[m2&(np.abs(dm2)<np.abs(dm1))]= ii[m2&(np.abs(dm2)<np.abs(dm1)),1]
    mask[m3&(np.abs(dm3)<np.abs(dm1))&(np.abs(dm3)<np.abs(dm2))]= ii[m3&(np.abs(dm3)<np.abs(dm1))&(np.abs(dm3)<np.abs(dm2)),2]
    return m[np.unique(mask[mask>=0].astype(int))], t[mask>=0]

def read_det_catalog(glob_query):
    start  = 0
    det = None
    for i,f in enumerate(np.sort(glob.glob(glob_query))):
        try:
            tmp = fio.FITS(f)[-1].read(columns=['mag_auto_J129','mag_auto_F184','mag_auto_H158','mag_auto_Y106', 'number', 'alphawin_j2000', 'deltawin_j2000','flux_auto','fluxerr_auto','flags'])
            # tmp = fio.FITS(f)[-1].read()
            if det is None:
                print(tmp.dtype)
                det = np.zeros(100000000,dtype=tmp.dtype)
            for col in det.dtype.names:
                det[col][start:start+len(tmp)] = tmp[col]
            start+=len(tmp)
        except:
            # print('-----fail'+f)
            pass
    return det


def read_truth_catalog(glob_query):
    dc2_coaddlist = fio.FITS('/hpc/group/cosmology/phy-lsst/public/dc2_sim_output/truth/coadd/dc2_coaddlist.fits.gz')[-1].read()
    dc2_coaddlist_dict = {row[0]: row for row in dc2_coaddlist}
    start  = 0
    truth = None
    for i,f in enumerate(np.sort(glob.glob(glob_query))):
        # print(i)
        try:
            tmp = fio.FITS(f)[-1].read(columns=['ra','dec', 'mag_J129','mag_F184','mag_H158','mag_Y106', 'ind', 'gal_star','x'])
            # tmp = fio.FITS(f)[-1].read()
            if truth is None:
                print(tmp.dtype)
                truth = np.zeros(100000000,dtype=tmp.dtype)
            # unique = extract_unique_entries(tmp)
            # del(tmp)
            mask = dc2_coaddlist_dict[extract_coordinates_regex(f)]
            # print(max(tmp['ra']), min(tmp['ra']))
            tmp = tmp[np.logical_and(np.logical_and(tmp['ra'] < mask['coadd_ra']+mask['d_ra'], tmp['ra'] > mask['coadd_ra']-mask['d_ra']),
                                     np.logical_and(tmp['dec'] < mask['coadd_dec']+mask['d_dec'], tmp['dec'] > mask['coadd_dec']-mask['d_dec']))]


            # print(max(unique['ra']), min(unique['ra']))
            # unique = unique[np.logical_and(np.logical_and(unique['ra'] < mask['coadd_ra']+mask['d_ra']/2, unique['ra'] > mask['coadd_ra']-mask['d_ra']/2),
            #                          np.logical_and(unique['dec'] < mask['coadd_dec']+mask['d_dec']/2, unique['dec'] > mask['coadd_dec']-mask['d_dec']/2))]
            
            for col in truth.dtype.names:
                truth[col][start:start+len(tmp)] = tmp[col]
            start+=len(tmp)
        except:
            # print('-----fail'+f)
            pass
    return truth

def read_truth_catalog_old(glob_query):
    start  = 0
    truth = None
    for i,f in enumerate(np.sort(glob.glob(glob_query))):
        # print(i)
        try:
            tmp = fio.FITS(f)[-1].read(columns=['ra','dec', 'mag_J129','mag_F184','mag_H158','mag_Y106', 'ind', 'gal_star','x'])
            # tmp = fio.FITS(f)[-1].read()
            if truth is None:
                print(tmp.dtype)
                truth = np.zeros(100000000,dtype=tmp.dtype)

            
            for col in truth.dtype.names:
                truth[col][start:start+len(tmp)] = tmp[col]
            start+=len(tmp)
        except:
            # print('-----fail'+f)
            pass
    return truth

# find redshift distribution for dc2
def nz_dist(file_path, ind):
    dc2_truth_gal = fio.FITS(file_path)[-1].read() # read in truth data for galaxy
    return dc2_truth_gal[ind]['z']

# # calculate threshold

# threshold = 0.85
# mag_threshold = np.zeros(4)
# # Plot the histogram
# for i in range(4):
# 	n1, bins1, patches1 = plt.hist(rtd['mag_'+fr[i]], bins=np.linspace(10,35,1300), histtype='step', label='All Input')
# 	n2, bins2, patches2 = plt.hist(rd['mag_auto_'+fr[i]], bins=np.linspace(10,35,1300), histtype='step', label='Detected')
# 	for j in range(n1.size):
# 		if(n2[j] > n1[j]*threshold and n2[j] > 1000):
# 			mag_threshold[i] = bins1[j]
# 	print('The threshold at ' + fr[i] + ' is ', mag_threshold[i])
# 	plt.clf()

def get_cov(filepath):

	data = np.genfromtxt(filepath)
	ndata = int(np.max(data[:,0]))+1

	print("Dimension of cov: %dx%d"%(ndata,ndata))

	ndata_min = int(np.min(data[:,0]))
	cov_g = np.zeros((ndata,ndata))
	cov_ng = np.zeros((ndata,ndata))
	for i in range(0,data.shape[0]):
		cov_g[int(data[i,0]),int(data[i,1])] =data[i,8]
		cov_g[int(data[i,1]),int(data[i,0])] =data[i,8]
		cov_ng[int(data[i,0]),int(data[i,1])] =data[i,9]
		cov_ng[int(data[i,1]),int(data[i,0])] =data[i,9]

	return cov_g, cov_ng, ndata

def load_cov(path):
    import os
    try:
        covdata = np.loadtxt(path)
    except:
        print('Skipping covariance, since output file missing.')
        return

    # Replace theta values with bin numbers.
    theta=np.sort(np.unique(covdata[:,2]))
    for i in range(len(theta)):
        covdata[np.where(covdata[:,2]==theta[i])[0],2]=i
        covdata[np.where(covdata[:,3]==theta[i])[0],3]=i

    # Populate covariance matrix.
    ndata=int(np.max(covdata[:,0]))+1
    ndata2=int(np.max(covdata[:,1]))+1
    assert ndata==ndata2

    cov=np.zeros((ndata,ndata))
    cov[:,:] = 0.0
    for i in range(0,covdata.shape[0]):
        cov[int(covdata[i,0]),int(covdata[i,1])]=covdata[i,8]
        cov[int(covdata[i,1]),int(covdata[i,0])]=covdata[i,8]
    covNG=np.zeros((ndata,ndata))
    covNG[:,:] = 0.0
    for i in range(0,covdata.shape[0]):
        covNG[int(covdata[i,0]),int(covdata[i,1])]=covdata[i,8]+covdata[i,9]
        covNG[int(covdata[i,1]),int(covdata[i,0])]=covdata[i,8]+covdata[i,9]

    covmat = (cov,covNG)
    return covmat