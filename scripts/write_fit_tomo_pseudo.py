from __future__ import print_function, division
import numpy as np
import twopoint
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

obj_type = 'roman'
cov_type = 'roman'
lens_type = 'len_'+obj_type
source_type = 'source_'+obj_type
bin_num = 5

gg_file_path = '/hpc/group/cosmology/denniswu/data/gg_corr_all.fits'
ng_file_path = '/hpc/group/cosmology/denniswu/data/ng_corr_21.fits'
ng_rand_file_path = '/hpc/group/cosmology/denniswu/data/ng_rand_21.fits'
nn_file_path = '/hpc/group/cosmology/denniswu/data/nn_corr_21.fits'
cov_file_path = '/hpc/group/cosmology/denniswu/covmat/'+obj_type+'/cov_'+obj_type
full_file_path = '/hpc/group/cosmology/denniswu/data/pseudo_'+obj_type+'.fits'
full_cov_file_path = '/hpc/group/cosmology/denniswu/cosmosis_files/data/pseudo_'+obj_type+'.fits'
# cov_file_path = '/hpc/group/cosmology/denniswu/covmat/'+cov_type+'/cov_'+cov_type
# full_file_path = '/hpc/group/cosmology/denniswu/data/pseudo_'+obj_type+'_'+cov_type+'.fits'
# full_cov_file_path = '/hpc/group/cosmology/denniswu/cosmosis_files/data/pseudo_'+obj_type+'_'+cov_type+'.fits'

bins_tomo = [1,2,3,4,5]
nbins_tomo = len(bins_tomo)

nbins_nz = 149

nbins = 20
min_sep, max_sep = 1, 100
bin_size = math.log(max_sep / min_sep)/ nbins        
logr = np.linspace(0, nbins*bin_size, nbins, endpoint=False, dtype=float)
logr += math.log(min_sep) + 0.5*bin_size
    
rnom = np.exp(logr)
half_bin = np.exp(0.5*bin_size)
left_edges = rnom / half_bin
right_edges = rnom * half_bin

print(left_edges)
print(right_edges)


''' Write lens and source nofz to twopoint fits file. '''
# Fetch nz dist, create hist
nz_dir = '/hpc/group/cosmology/denniswu/hist_file.npy'
nz_source = np.load(nz_dir,allow_pickle='TRUE').item()[source_type]
nz_lens = np.load(nz_dir,allow_pickle='TRUE').item()[lens_type]

hist_source, bins_source = zip(*nz_source)
hist_lens, bins_lens = zip(*nz_lens)
bins_source = bins_source[0]
bins_lens = bins_lens[0]

n_source = [np.sum(arr) for arr in list(hist_source)]
n_lens = [np.sum(arr) for arr in list(hist_lens)]

bins_source_mid = [(bins_source[i]+bins_source[i+1])/2 for i in range(nbins_nz)]
bins_lens_mid = [(bins_lens[i]+bins_lens[i+1])/2 for i in range(nbins_nz)]

# input nofz info
nz_source = twopoint.NumberDensity(
    NOFZ_NAMES[0],
    bins_source[:-1],
    bins_source_mid,
    bins_source[1:],
    hist_source)

nz_source.ngal      = n_source
nz_source.sigma_e   = np.repeat([0.379411],nbins_tomo)
nz_source.area      = 20

nz_lens = twopoint.NumberDensity(
    NOFZ_NAMES[1],
    bins_lens[:-1],
    bins_lens_mid,
    bins_lens[1:],
    hist_lens)

nz_lens.ngal = n_lens
nz_lens.sigma_e   = np.repeat([0.379411],nbins_tomo)
nz_lens.area = 20

kernels  = [nz_source, nz_lens]


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

exts[0].angular_bin = np.tile(np.arange(length),np.sum(bins_tomo))
exts[0].angle_min   = np.tile(left_edges, np.sum(bins_tomo))
exts[0].angle_max   = np.tile(right_edges, np.sum(bins_tomo))
exts[0].angle       = np.tile(gg_data.meanr, np.sum(bins_tomo))         
exts[0].bin1        = np.repeat(bins_tomo, [x * nbins for x in reversed(range(1, nbins_tomo + 1))])
exts[0].bin2        = np.repeat(np.concatenate([bins_tomo[i:] for i in range(nbins_tomo)]), nbins)
exts[0].value       = np.tile(gg_data.xip, np.sum(bins_tomo))
exts[0].npairs      = np.tile(gg_data.npairs, np.sum(bins_tomo))
exts[0].weight      = np.tile(gg_data.weight, np.sum(bins_tomo))

exts[1].angular_bin = np.tile(np.arange(length),np.sum(bins_tomo))
exts[1].angle_min   = np.tile(left_edges, np.sum(bins_tomo))
exts[1].angle_max   = np.tile(right_edges, np.sum(bins_tomo))
exts[1].angle       = np.tile(gg_data.meanr, np.sum(bins_tomo))           
exts[1].bin1        = np.repeat(bins_tomo, [x * nbins for x in reversed(range(1, nbins_tomo + 1))])
exts[1].bin2        = np.repeat(np.concatenate([bins_tomo[i:] for i in range(nbins_tomo)]), nbins)
exts[1].value       = np.tile(gg_data.xim, np.sum(bins_tomo))
exts[1].npairs      = np.tile(gg_data.npairs, np.sum(bins_tomo))
exts[1].weight      = np.tile(gg_data.weight, np.sum(bins_tomo))

# tangential shear, shear-pos (+, -)
ng_data = read_treecorr(ng_file_path, 'ng')
ng_rand_data = read_treecorr(ng_rand_file_path, 'ng')
length = len(ng_data.meanr)

gammat_size = {'truth': 19, 'rubin': 24, 'roman': 24}
gammat_bin1 = {'truth': [1,1,1,1,1,2,2,2,2,2,3,3,3,3,4,4,4,5,5], 'rubin': [1,1,1,1,1,2,2,2,2,2,3,3,3,3,3,4,4,4,4,4,5,5,5,5], 'roman': [1,1,1,1,1,2,2,2,2,2,3,3,3,3,3,4,4,4,4,4,5,5,5,5]}
gammat_bin2 = {'truth': [1,2,3,4,5,1,2,3,4,5,2,3,4,5,3,4,5,4,5], 'rubin': [1,2,3,4,5,1,2,3,4,5,1,2,3,4,5,1,2,3,4,5,2,3,4,5], 'roman': [1,2,3,4,5,1,2,3,4,5,1,2,3,4,5,1,2,3,4,5,2,3,4,5]}

exts[2].angular_bin   = np.tile(np.arange(length), gammat_size[obj_type])
exts[2].angle_min     = np.tile(left_edges, gammat_size[obj_type])
exts[2].angle_max     = np.tile(right_edges, gammat_size[obj_type])                       
exts[2].angle         = np.tile(ng_data.meanr, gammat_size[obj_type])        
exts[2].bin1          = np.repeat(gammat_bin1[obj_type], nbins)
exts[2].bin2          = np.repeat(gammat_bin2[obj_type], nbins)
exts[2].value         = np.tile(ng_data.xi, gammat_size[obj_type])
exts[2].npairs        = np.tile(ng_data.npairs, gammat_size[obj_type])
exts[2].weight        = np.tile(ng_data.weight, gammat_size[obj_type])
exts[2].random_npairs = np.tile(ng_rand_data.npairs, gammat_size[obj_type])
exts[2].random_weight = np.tile(ng_rand_data.weight, gammat_size[obj_type])


exts[3].angular_bin   = np.tile(np.arange(length),gammat_size[obj_type])
exts[3].angle_min     = np.tile(left_edges, gammat_size[obj_type])
exts[3].angle_max     = np.tile(right_edges, gammat_size[obj_type])       
exts[3].angle         = np.tile(ng_data.meanr, gammat_size[obj_type])        
exts[3].bin1          = np.repeat(gammat_bin1[obj_type], nbins)
exts[3].bin2          = np.repeat(gammat_bin2[obj_type], nbins)
exts[3].value         = np.tile(ng_data.xi_im, gammat_size[obj_type])
exts[3].npairs        = np.tile(ng_data.npairs, gammat_size[obj_type])
exts[3].weight        = np.tile(ng_data.weight, gammat_size[obj_type])
exts[3].random_npairs = np.tile(ng_rand_data.npairs, gammat_size[obj_type])
exts[3].random_weight = np.tile(ng_rand_data.weight, gammat_size[obj_type])

# galaxy clustering, pos-pos
nn_data = read_treecorr(nn_file_path, 'nn')
length = len(nn_data.meanr)

exts[4].angular_bin   = np.tile(np.arange(length),nbins_tomo)
exts[4].angle_min     = np.tile(left_edges, nbins_tomo)
exts[4].angle_max     = np.tile(right_edges, nbins_tomo)       
exts[4].angle         = np.tile(nn_data.meanr, nbins_tomo)           
exts[4].bin1          = np.repeat(bins_tomo, nbins)
exts[4].bin2          = np.repeat(bins_tomo, nbins)
exts[4].value         = np.tile(nn_data.xi, nbins_tomo)
exts[4].npairs        = np.tile(nn_data.npairs, nbins_tomo)
exts[4].weight        = np.tile(nn_data.weight, nbins_tomo)
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

covmat = util.load_cov(cov_file_path)	

if covmat is not None:
    lengths = [len(s) for s in fits.spectra]
    print(lengths)
    names = [s.name for s in fits.spectra]
    print(names)

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

T1 = twopoint.TwoPointFile.from_fits(full_cov_file_path)
T1.plots(full_cov_file_path[:-5], latex = False)

# length=get_cov_lengths(fits)

# fits.covmat_info=twopoint.CovarianceMatrixInfo('COVMAT',TWO_POINT_NAMES,length,gg_data.cov)
# fits.to_fits(g_cov_file_path,clobber=True)
# fits.covmat_info=twopoint.CovarianceMatrixInfo('COVMAT',TWO_POINT_NAMES,length,ng_data.cov)
# fits.to_fits(ng_cov_file_path,clobber=True)


