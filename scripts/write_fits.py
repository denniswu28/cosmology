from __future__ import print_function, division
import scripts.twopoint as twopoint
import numpy as np
import treecorr
import utilities as util
import math
from astropy.io import fits as fitsio

TWO_POINT_NAMES = ['xip','xim','gammat','gammax','wtheta']
TWO_POINT_NAMES_MOD = ['xip','xim','gammat','wtheta']
NOFZ_NAMES = ['nz_source', 'nz_lens']
COV_NAME = 'COVMAT'
fr = ['Y106','J129','H158','F184']

gg_file_path = '/hpc/group/cosmology/denniswu/data/gg_corr_all.fits'
ng_file_path = '/hpc/group/cosmology/denniswu/data/ng_corr_21.fits'
ng_rand_file_path = '/hpc/group/cosmology/denniswu/data/ng_rand_21.fits'
nn_file_path = '/hpc/group/cosmology/denniswu/data/nn_corr_21.fits'
g_cov_file_path = '/hpc/group/cosmology/denniswu/covmat/cov_dc2_21'
full_file_path = '/hpc/group/cosmology/denniswu/data/2pt_CosmoDC2_21_det_shear.fits'
full_cov_file_path = '/hpc/group/cosmology/denniswu/cosmosis_files/data/2pt_all_CosmoDC2_21_det_shear.fits'

mag_threshold = 21
dec_min = -42
dec_max = -38
ra_min = 51
ra_max = 56

nbins = 20
min_sep, max_sep = 1, 100
bin_size = math.log(max_sep / min_sep)/ nbins        
logr = np.linspace(0, nbins*bin_size, nbins, endpoint=False, dtype=float)
logr += math.log(min_sep) + 0.5*bin_size
    
rnom = np.exp(logr)
half_bin = np.exp(0.5*bin_size)
left_edges = rnom / half_bin
right_edges = rnom * half_bin


''' Write a dummy lens and source nofz to twopoint fits file. '''
# Read in and cut dectection objects
dc2_det = util.read_det_catalog('/hpc/group/cosmology/phy-lsst/public/dc2_sim_output/detection/dc2_det_*.fits.gz')
dc2_det = dc2_det[dc2_det['number']>0] # select valid
dc2_det = dc2_det[(dc2_det['flags']<4)&(dc2_det['flux_auto']/dc2_det['fluxerr_auto']>18)] # select flag=0 and flux/err ratio >5
dc2_det = dc2_det[(dc2_det['alphawin_j2000']>ra_min)&(dc2_det['alphawin_j2000']<ra_max)&(dc2_det['deltawin_j2000']>dec_min)&(dc2_det['deltawin_j2000']<dec_max)] # boundary limit

# Read in and cut truth objects
dc2_truth = util.read_truth_catalog('/hpc/group/cosmology/phy-lsst/public/dc2_sim_output/truth/coadd/dc2_index_*.fits.gz')
dc2_truth = dc2_truth[dc2_truth['x']>0] # select valid
dc2_truth = dc2_truth[dc2_truth['gal_star'] == 0] # select galaxy
dc2_truth = dc2_truth[np.logical_and(np.logical_and(dc2_truth['ra'] < ra_max, dc2_truth['ra'] > ra_min),np.logical_and(dc2_truth['dec'] < dec_max, dc2_truth['dec'] > dec_min))]

# Match both catalogs
dc2_det_match, dc2_truth_match = util.get_match(dc2_det, dc2_truth)

# Single out the lens
dc2_truth_thre = dc2_truth_match[dc2_truth_match['mag_'+fr[2]] < mag_threshold]

# Output the source, lens number
print("Matched source number: ",dc2_truth_match['ra'].size)
print("Matched lens number w/ m_H<" + str(mag_threshold) + " : ",dc2_truth_thre['ra'].size)

# Fetch nz dist from truth, create hist
ddir = '/hpc/group/cosmology/phy-lsst/public/dc2_sim_output/truth/dc2_truth_gal.fits'
nz_lens = util.nz_dist(ddir, dc2_truth_thre['ind'])
nz_source = util.nz_dist(ddir, dc2_truth_match['ind'])

hist_lens, bins_lens = np.histogram(nz_lens, bins = 300, range=(0,3), density = True)
bins_lens_mid = [(bins_lens[i]+bins_lens[i+1])/2 for i in range(300)]
print(bins_lens_mid[0:10]) 
hist_source, bins_source = np.histogram(nz_source, bins = 300, range=(0,3), density = True)
bins_source_mid = [(bins_source[i]+bins_source[i+1])/2 for i in range(300)]
print(bins_source_mid[0:10]) 

# input nofz info
nz_source = twopoint.NumberDensity(
    NOFZ_NAMES[0],
    bins_source[:-1],
    bins_source_mid,
    bins_source[1:],
    [hist_source])

nz_source.ngal      = [dc2_truth_match['ra'].size]
nz_source.sigma_e   = [0.379411]
nz_source.area      = 20

nz_lens = twopoint.NumberDensity(
    NOFZ_NAMES[1],
    bins_lens[:-1],
    bins_lens_mid,
    bins_lens[1:],
    [hist_lens])

nz_lens.ngal = [dc2_truth_thre['ra'].size]
nz_lens.area = 20
kernels  = (nz_source, nz_lens)


''' Set up the 2pt class '''
dtype1=[twopoint.Types.galaxy_shear_plus_real,twopoint.Types.galaxy_shear_minus_real,twopoint.Types.galaxy_position_real,twopoint.Types.galaxy_position_real,twopoint.Types.galaxy_position_real]
dtype2=[twopoint.Types.galaxy_shear_plus_real,twopoint.Types.galaxy_shear_minus_real,twopoint.Types.galaxy_shear_plus_real,twopoint.Types.galaxy_position_real,twopoint.Types.galaxy_position_real]
nznameindex1=[0,0,1,1,1]
nznameindex2=[0,0,0,0,1]

def get_length(n,n2=None):
    if n2 is not None:
        return n*n2
    else:
        return n*(n+1)/2

# Read 2pt data from default treecorr fits file
def read_treecorr(file_path, corr_type):
    if corr_type == 'gg':
        data = treecorr.GGCorrelation(min_sep=min_sep, max_sep=max_sep, nbins=nbins, sep_units='arcmin', var_method='jackknife')
        data.read(file_path, file_type="FITS")
    if corr_type == 'ng':
        data = treecorr.NGCorrelation(min_sep=min_sep, max_sep=max_sep, nbins=nbins, sep_units='arcmin', var_method='jackknife')
        data.read(file_path, file_type="FITS")
    if corr_type == 'nn':
        data = treecorr.NNCorrelation(min_sep=min_sep, max_sep=max_sep, nbins=nbins, sep_units='arcmin', var_method='jackknife')
        data.read(file_path, file_type="FITS")
    return data

# def get_cov_lengths(covpath):
#     # Calculate length of covariance blocks. Exception for gammat, which reads a file that stores which bin combinations have been rejected in the covariance calcultion.
#     # Make cov lengths
#     _, _, ndata = util.get_cov(covpath)
#     return ndata
    
    
exts = []

for i,name in enumerate(TWO_POINT_NAMES):
    exts.append(twopoint.SpectrumMeasurement(
        name, # hdu name
        ([],[]), # tomographic bins
        (dtype1[i], dtype2[i]), # type of 2pt statistic
        (NOFZ_NAMES[nznameindex1[i]], NOFZ_NAMES[nznameindex2[i]]), # associated nofz
        "SAMPLE", # window function
        None, # id
        None, # value
        npairs=None, # pair counts
        angle=None,
        angle_unit='arcmin')) # units

# cosmic shear, shear-shear
gg_data = read_treecorr(gg_file_path, 'gg')
length = len(gg_data.meanr)

exts[0].angular_bin = np.arange(length)
exts[0].angle_min   = left_edges
exts[0].angle_max   = right_edges
exts[0].angle       = gg_data.meanr           
exts[0].bin1        = np.ones(length)
exts[0].bin2        = np.ones(length)
exts[0].value       = gg_data.xip
exts[0].npairs      = gg_data.npairs
exts[0].weight      = gg_data.weight
exts[1].angular_bin = np.arange(length)
exts[1].angle_min   = left_edges
exts[1].angle_max   = right_edges
exts[1].angle       = np.exp(gg_data.meanlogr)              
exts[1].bin1        = np.ones(length)
exts[1].bin2        = np.ones(length)
exts[1].value       = gg_data.xim
exts[1].npairs      = gg_data.npairs
exts[1].weight      = gg_data.weight


# tangential shear, shear-pos
ng_data = read_treecorr(ng_file_path, 'ng')
ng_rand_data = read_treecorr(ng_rand_file_path, 'ng')
length = len(ng_data.meanr)

exts[2].angular_bin   = np.arange(length)
exts[2].angle_min     = left_edges
exts[2].angle_max     = right_edges                        
exts[2].angle         = ng_data.meanr
exts[2].bin1          = np.ones(length)
exts[2].bin2          = np.ones(length)
exts[2].value         = ng_data.xi
exts[2].npairs        = ng_data.npairs
exts[2].weight        = ng_data.weight
exts[2].random_npairs = ng_rand_data.npairs
exts[2].random_weight = ng_rand_data.weight

exts[3].angular_bin   = np.arange(length)
exts[3].angle_min     = left_edges
exts[3].angle_max     = right_edges                        
exts[3].angle         = ng_data.meanr
exts[3].bin1          = np.ones(length)
exts[3].bin2          = np.ones(length)
exts[3].value         = ng_data.xi_im
exts[3].npairs        = ng_data.npairs
exts[3].weight        = ng_data.weight
exts[3].random_npairs = ng_rand_data.npairs
exts[3].random_weight = ng_rand_data.weight

# galaxy clustering, pos-pos
nn_data = read_treecorr(nn_file_path, 'nn')
length = len(nn_data.meanr)

exts[4].angular_bin   = np.arange(length)
exts[4].angle_min     = left_edges
exts[4].angle_max     = right_edges            
exts[4].angle         = nn_data.meanr           
exts[4].bin1          = np.ones(length)
exts[4].bin2          = np.ones(length)
exts[4].value         = nn_data.xi
exts[4].npairs        = nn_data.npairs
exts[4].weight        = nn_data.weight
# exts[4].random_npairs = np.zeros(length)
# exts[4].random_weight = np.zeros(length)
# exts[4].dr_npairs     = np.zeros(length)
# exts[4].dr_weight     = np.zeros(length)
# exts[4].rd_npairs     = np.zeros(length)
# exts[4].rd_weight     = np.zeros(length)


''' Write file without covariance (all data vectors) '''
fits = twopoint.TwoPointFile(exts, kernels, None, None)
# fits.to_fits(full_file_path, clobber=True)

# print(TWO_POINT_NAMES)
# print(fits.spectra[3]) #[<Spectrum: xip>, <Spectrum: xim>, <Spectrum: gammat>, <Spectrum: gammax>, <Spectrum: wtheta>]

del fits.spectra[3]

print('spectra:', fits.spectra)

# print(fits.spectra[-1])
# print(fits.spectra[-1].name)

'''Writes the covariance info into a covariance object and saves to new fits files.'''

covfile = '/hpc/group/cosmology/denniswu/covmat/cov_dc2_7'
covmat = util.load_cov(covfile)	

if covmat is not None:
    lengths = [len(s) for s in fits.spectra]
    names = [s.name for s in fits.spectra]

# Writes the covariance info into a covariance object and saves to 2point fits file.
fits.covmat_info=twopoint.CovarianceMatrixInfo('COVMAT', TWO_POINT_NAMES_MOD, lengths, covmat[0])
fits.to_fits(full_cov_file_path, clobber=True)


fitsfile = fitsio.open(full_cov_file_path)
extension = fitsfile['COVMAT']
covmat_info = twopoint.CovarianceMatrixInfo.from_fits(extension)
print('covmat_info names:', covmat_info.names)
for extension in fitsfile:
    if extension.header.get('2PTDATA'):
        name = extension.name
        print('extension name:',name)
        






# length=get_cov_lengths(fits)

# fits.covmat_info=twopoint.CovarianceMatrixInfo('COVMAT',TWO_POINT_NAMES,length,gg_data.cov)
# fits.to_fits(g_cov_file_path,clobber=True)
# fits.covmat_info=twopoint.CovarianceMatrixInfo('COVMAT',TWO_POINT_NAMES,length,ng_data.cov)
# fits.to_fits(ng_cov_file_path,clobber=True)


