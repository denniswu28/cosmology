import numpy as np
from numpy import linalg as LA
import matplotlib.pyplot as plt
import sys

def get_cov(filename):

	data = np.genfromtxt(filename)
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


covfile = '/hpc/group/cosmology/denniswu/covmat/cov_dc2_21'


'''Plot the covariance matrix'''
c_g, c_ng, ndata = get_cov(covfile)	

print(c_g[0:10,0:10])

cov = c_g

b = np.sort(LA.eigvals(cov))
print("min+max eigenvalues cov: %e, %e"%(np.min(b), np.max(b)))
if(np.min(b)<=0.):
	print("non-positive eigenvalue encountered! Covariance Invalid!")
	exit()

print("Covariance is postive definite!")

cmap = 'seismic'

pp_norm = np.zeros((ndata,ndata))
for i in range(ndata):
	for j in range(ndata):
		pp_norm[i][j] = cov[i][j]/ np.sqrt(cov[i][i]*cov[j][j])

print(pp_norm[0:10,0:10])

print("Plotting correlation matrix ...")

plot_path = covfile+'_plot.pdf'

fig = plt.figure()
ax = fig.add_subplot(1, 1, 1)
plt.axvline(x=20,color='black')
plt.axvline(x=40,color='black')
plt.axvline(x=60,color='black')
plt.axhline(y=20,color='black')
plt.axhline(y=40,color='black')
plt.axhline(y=60,color='black')
ax.invert_yaxis()
im3 = ax.imshow(pp_norm, cmap=cmap, vmin=-1, vmax=1)
fig.colorbar(im3, orientation='vertical')

# ax.text(15, 15, , fontsize=14)
ax.annotate(r'$\xi_+^{ij}(\theta)$', xy=(15.0, -15.0), xycoords='data', annotation_clip=False)

# ax.text(265, -15, r'$\xi_-^{ij}(\theta)$', fontsize=14)
# ax.text(565, -15, r'$\gamma_t^{ij}(\theta)$', fontsize=14)
# ax.text(815, -15, r'$w^{i}(\theta)$', fontsize=14)

# ax.text(905, 95, r'$\xi_+$', fontsize=14)
# ax.text(905, 295, r'$\xi_-$', fontsize=14)
# ax.text(905, 595, r'$\gamma_t$', fontsize=14)
# ax.text(905, 845, r'$w$', fontsize=14)

plt.savefig(plot_path,dpi=2000)
plt.show()
print("Plot saved as %s"%(plot_path))