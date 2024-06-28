from __future__ import print_function, division
import numpy as np
import matplotlib
import twopoint
import numpy as np
import pandas as pd
import treecorr
import utilities as util
import math
from astropy.io import fits as fitsio

matplotlib.rcParams.update(matplotlib.rcParamsDefault)

TWO_POINT_NAMES = ['xip','xim','gammat','gammax','wtheta']
TWO_POINT_NAMES_MOD = ['xip','xim','gammat','wtheta']
NOFZ_NAMES = ['nz_source', 'nz_lens']
COV_NAME = 'COVMAT'
fr = ['Y106','J129','H158','F184']

obj_type = 'roman'
cov_type = 'roman'
nz_type = 'roman'
lens_type = 'len_'+nz_type
source_type = 'source_'+nz_type
bin_num = 5
vector_dir_path = '/hpc/group/cosmology/denniswu/cosmosis_files/output/'+obj_type+'_'+obj_type+'/data_vector/'
cov_file_path = '/hpc/group/cosmology/denniswu/covmat/'+cov_type+'/cov_'+cov_type
full_file_path = '/hpc/group/cosmology/denniswu/data/'+obj_type+'_'+cov_type+'.fits'
full_cov_file_path = '/hpc/group/cosmology/denniswu/cosmosis_files/data/'+obj_type+'_'+cov_type+'.fits'

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


''' Write lens and source nofz to twopoint fits file. '''
# Fetch nz dist from truth, create hist
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

# Read 2pt data from data_vector text files
data_vectors = pd.DataFrame()
vector_filenames = ['2pt_angle.txt','2pt_theory.txt','2pt_bin1.txt','2pt_bin2.txt']
for filename in vector_filenames:
    data = np.loadtxt(vector_dir_path+filename, comments = "#")
    column_label = filename.replace('.txt', '')  # Strip the file extension for the label
    data_vectors[column_label] = data
twopt_indices = {'truth': {
    'xip': 0,
    'xim': 300,
    'gammat': 600,
    'wtheta': 980},
    'rubin': {
    'xip': 0,
    'xim': 300,
    'gammat': 600,
    'wtheta': 1080},
    'roman': {
    'xip': 0,
    'xim': 300,
    'gammat': 600,
    'wtheta': 1080},
}

def assign_type(index):
    # Sort the dictionary by value to ensure the indices are in order
    sorted_types = sorted(twopt_indices[obj_type].items(), key=lambda x: x[1])
    for i in range(len(sorted_types)):
        # Check if the current index is less than the starting index of the next type
        if i + 1 < len(sorted_types) and index < sorted_types[i + 1][1]:
            return sorted_types[i][0]
        # For the last type segment
        elif i + 1 == len(sorted_types):
            return sorted_types[i][0]
    return "Unknown"
data_vectors['2pt_name'] = [assign_type(i) for i in range(len(data_vectors))]
print(data_vectors.iloc[0:5])
print(data_vectors.iloc[980:985])

# Setting up the 2pt classes and loading the 2pt info
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
xip_data = data_vectors[data_vectors['2pt_name'] == 'xip']
xim_data = data_vectors[data_vectors['2pt_name'] == 'xim']
length = nbins

exts[0].angular_bin = np.tile(np.arange(length),np.sum(bins_tomo))
exts[0].angle_min   = np.tile(left_edges, np.sum(bins_tomo))
exts[0].angle_max   = np.tile(right_edges, np.sum(bins_tomo))
exts[0].angle       = xip_data['2pt_angle']         
exts[0].bin1        = xip_data['2pt_bin1']
exts[0].bin2        = xip_data['2pt_bin2']
exts[0].value       = xip_data['2pt_theory']
# exts[0].npairs      = np.tile(gg_data.npairs, np.sum(bins_tomo))
# exts[0].weight      = np.tile(gg_data.weight, np.sum(bins_tomo))

exts[1].angular_bin = np.tile(np.arange(length),np.sum(bins_tomo))
exts[1].angle_min   = np.tile(left_edges, np.sum(bins_tomo))
exts[1].angle_max   = np.tile(right_edges, np.sum(bins_tomo))
exts[1].angle       = xim_data['2pt_angle']         
exts[1].bin1        = xim_data['2pt_bin1']
exts[1].bin2        = xim_data['2pt_bin2']
exts[1].value       = xim_data['2pt_theory']
# exts[1].npairs      = np.tile(gg_data.npairs, np.sum(bins_tomo))
# exts[1].weight      = np.tile(gg_data.weight, np.sum(bins_tomo))

# tangential shear, shear-pos (+, -)
gammat_data = data_vectors[data_vectors['2pt_name'] == 'gammat']
length = nbins

gammat_size = {'truth': 19, 'rubin': 24, 'roman': 24}
gammat_bin1 = {'truth': [1,1,1,1,1,2,2,2,2,2,3,3,3,3,4,4,4,5,5], 'rubin': [1,1,1,1,1,2,2,2,2,2,3,3,3,3,3,4,4,4,4,4,5,5,5,5], 'roman': [1,1,1,1,1,2,2,2,2,2,3,3,3,3,3,4,4,4,4,4,5,5,5,5]}
gammat_bin2 = {'truth': [1,2,3,4,5,1,2,3,4,5,2,3,4,5,3,4,5,4,5], 'rubin': [1,2,3,4,5,1,2,3,4,5,1,2,3,4,5,1,2,3,4,5,2,3,4,5], 'roman': [1,2,3,4,5,1,2,3,4,5,1,2,3,4,5,1,2,3,4,5,2,3,4,5]}

exts[2].angular_bin   = np.tile(np.arange(length),gammat_size[cov_type])
exts[2].angle_min     = np.tile(left_edges, gammat_size[cov_type])
exts[2].angle_max     = np.tile(right_edges, gammat_size[cov_type])                       
exts[2].angle         = gammat_data['2pt_angle'][:gammat_size[cov_type]]
exts[2].bin1          = gammat_data['2pt_bin1'][:gammat_size[cov_type]]
exts[2].bin2          = gammat_data['2pt_bin2'][:gammat_size[cov_type]]
exts[2].value         = gammat_data['2pt_theory'][:gammat_size[cov_type]]
# exts[2].npairs        = np.tile(ng_data.npairs, 19)
# exts[2].weight        = np.tile(ng_data.weight, 19)
# exts[2].random_npairs = np.tile(ng_rand_data.npairs, 19)
# exts[2].random_weight = np.tile(ng_rand_data.weight, 19)


exts[3].angular_bin   = np.tile(np.arange(length),gammat_size[cov_type])
exts[3].angle_min     = np.tile(left_edges, gammat_size[cov_type])
exts[3].angle_max     = np.tile(right_edges, gammat_size[cov_type])                       
exts[3].angle         = gammat_data['2pt_angle'][:gammat_size[cov_type]]
exts[3].bin1          = gammat_data['2pt_bin1'][:gammat_size[cov_type]]
exts[3].bin2          = gammat_data['2pt_bin2'][:gammat_size[cov_type]]
exts[3].value         = gammat_data['2pt_theory'][:gammat_size[cov_type]]
# exts[3].npairs        = np.tile(ng_data.npairs, 19)
# exts[3].weight        = np.tile(ng_data.weight, 19)
# exts[3].random_npairs = np.tile(ng_rand_data.npairs, 19)
# exts[3].random_weight = np.tile(ng_rand_data.weight, 19)

# galaxy clustering, pos-pos
wtheta_data = data_vectors[data_vectors['2pt_name'] == 'wtheta']
length = nbins

exts[4].angular_bin   = np.tile(np.arange(length),nbins_tomo)
exts[4].angle_min     = np.tile(left_edges, nbins_tomo)
exts[4].angle_max     = np.tile(right_edges, nbins_tomo)       
exts[4].angle         = wtheta_data['2pt_angle']
exts[4].bin1          = wtheta_data['2pt_bin1']
exts[4].bin2          = wtheta_data['2pt_bin2']
exts[4].value         = wtheta_data['2pt_theory']
# exts[4].npairs        = wtheta_data['2pt_angle']
# exts[4].weight        = wtheta_data['2pt_angle']
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
print('spectra:', fits.spectra)

print(fits.spectra[3]) #[<Spectrum: xip>, <Spectrum: xim>, <Spectrum: gammat>, <Spectrum: gammax>, <Spectrum: wtheta>]

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

print('matplotlib: {}'.format(matplotlib.__version__))

T1 = twopoint.TwoPointFile.from_fits(full_cov_file_path)
T1.plots(full_cov_file_path[:-5], latex = False)

# length=get_cov_lengths(fits)

# fits.covmat_info=twopoint.CovarianceMatrixInfo('COVMAT',TWO_POINT_NAMES,length,gg_data.cov)
# fits.to_fits(g_cov_file_path,clobber=True)
# fits.covmat_info=twopoint.CovarianceMatrixInfo('COVMAT',TWO_POINT_NAMES,length,ng_data.cov)
# fits.to_fits(ng_cov_file_path,clobber=True)


