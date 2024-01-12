from __future__ import print_function, division
import numpy as np
import pandas as pd
import math
import matplotlib
matplotlib.use ('agg')
import matplotlib.pyplot as plt
from matplotlib.colors import LogNorm
# import seaborn as sns

TWO_POINT_NAMES = ['xip','xim','gammat','gammax','wtheta']
NOFZ_NAMES = ['nz_source', 'nz_lens']
COV_NAME = 'COVMAT'
fr = ['Y106','J129','H158','F184']

g_cov_file_path = '/hpc/group/cosmology/denniswu/covmat/dc2_cov'


'''Plot the covariance matrix'''

covdata = np.loadtxt(g_cov_file_path)
# cov_df = pd.DataFrame(covdata, columns = ['i','j', 'theta_1', 'theta_2', 'func_index', 'index_1', 'index_2', 'index_3', 'gaussian', 'n_gaussian'])

# Populate covariance matrix.
ndata=int(np.max(covdata[:,0]))+1
ndata2=int(np.max(covdata[:,1]))+1
assert ndata==ndata2

cov=np.zeros((ndata, ndata))
cov[:,:] = 0.0
for i in range(0, covdata.shape[0]):
    cov[int(covdata[i,0]),int(covdata[i,1])]=covdata[i,8]
    cov[int(covdata[i,1]),int(covdata[i,0])]=covdata[i,8]

cov = cov-np.min(cov)
# print(cov[1:20,1:20])
plt.imshow(cov, cmap='autumn', norm=LogNorm(vmin=np.min(cov), vmax=np.max(cov)))
plt.title('Covariance Matrix')
# plt.colorbar()
plt.savefig('corr_mat.pdf')
plt.clf()

# sns.heatmap(cov, annot=False, cmap=plt.cm.YlOrRd_r)
# plt.savefig('corr_mat_2.pdf')
# plt.clf()