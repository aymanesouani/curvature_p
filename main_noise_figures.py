import matplotlib
import trimesh
from sklearn.manifold import MDS

matplotlib.use('TkAgg')
import os
import generate_quadric_surfaces as sps
import slam.io as sio
import numpy as np
import matplotlib.pyplot as plt
from sklearn.metrics import mean_squared_error, median_absolute_error, mean_squared_log_error, r2_score
from scipy.stats import pearsonr

def corr(x, y):
	n = len(x)
	if n != len(y):
		raise ValueError('x and y must have the same length.')
	
	if n < 2:
		raise ValueError('x and y must have length at least 2.')
	
	x = np.asarray(x)
	y = np.asarray(y)
	# dtype is the data type for the calculations.  This expression ensures
	# that the data type is at least 64 bit floating point.  It might have
	# more precision if the input is, for example, np.longdouble.
	dtype = type(1.0 + x[0] + y[0])
	
	if n == 2:
		return dtype(np.sign(x[1] - x[0]) * np.sign(y[1] - y[0])), 1.0
	
	xmean = np.mean(x)
	ymean = np.mean(y)
	
	# By using `astype(dtype)`, we ensure that the intermediate calculations
	# use at least 64 bit floating point.
	xm = x - xmean
	ym = y - ymean
	
	# Unlike np.linalg.norm or the expression sqrt((xm*xm).sum()),
	# scipy.linalg.norm(xm) does not overflow if xm is, for example,
	# [-5e210, 5e210, 3e200, -3e200]
	normxm = np.linalg.norm(xm)
	normym = np.linalg.norm(ym)
	
	r = np.dot(np.transpose(xm / normxm), ym / normym)
	
	# Presumably, if abs(r) > 1, then it is only some small artifact of
	# floating point arithmetic.
	r = max(min(r, 1.0), -1.0)


	# As explained in the docstring, the p-value can be computed as
	#     p = 2*dist.cdf(-abs(r))
	# where dist is the beta distribution on [-1, 1] with shape parameters
	# a = b = n/2 - 1.  `special.btdtr` is the CDF for the beta distribution
	# on [0, 1].  To use it, we make the transformation  x = (r + 1)/2; the
	# shape parameters do not change.  Then -abs(r) used in `cdf(-abs(r))`
	# becomes x = (-abs(r) + 1)/2 = 0.5*(1 - abs(r)).  (r is cast to float64
	# to avoid a TypeError raised by btdtr when r is higher precision.)
	ab = n / 2 - 1
	return r
if __name__ == '__main__':
	curv_types = ['Peyre', 'Rusinkiewicz', 'PetitJean', 'Taubin', 'Dong']  # , 'Peyre', 'fem', 'boix', 'bary']
	curv_types_leg = curv_types.copy()
	curv_types_leg.append('analytic')
	meshes = list()
	couleur = ['red', 'green', 'yellow', 'purple', 'brown']
	Ks = [[0, 1], [1, 1], [-1, 1]]
	Nstep =[46, 32, 22, 19, 15]  # �r exemple
	Sigma = [0.2, 0.4, 0.5, 0.55, 0.6, 0.7, 0.8]
	list_exampels = list()
	os.chdir(r"/hpc/meca/users/souani.a/curvature/Compute_curvatures/Curvatures")
	
	for i in range(len(Ks)):
		for nstep in Nstep:
			fig = plt.figure()
			ax = fig.add_subplot(111)
			fig.suptitle("MSE function of Noise " +'hex_quad_k1_' + str(Ks[i][0]) + '_k2_' + str(Ks[i][1])+ '_nstep_' + str(nstep))
			plt.xlabel(' Noise')
			plt.ylabel(' MSE ')
			s = 0
			for curv_type in curv_types:
				err= list()
				for sigma in Sigma:
					cc = 'hex_quad_k1_' + str(Ks[i][0]) + '_k2_' + str(Ks[i][1]) + '_' + str(nstep) + '_noise_' + str(
						int(sigma * 1000))
					# os.chdir(r"/hpc/meca/users/souani.a/curvature/Compute_curvatures/Noise")
					# os.mkdir(cc )
					# os.chdir(r"/hpc/meca/users/souani.a/curvature/Compute_curvatures/Noise/"+cc)
					# sio.write_mesh(mesh, cc+'.gii')
					list_exampels.append(cc)
					mesh_file = 'hex_quad_k1_' + str(Ks[i][0]) + '_k2_' + str(Ks[i][1]) + '_' + str(nstep) + '.gii'
					os.chdir("/hpc/meca/users/souani.a/curvature/Compute_curvatures/Curvatures/0Noise")
					mesh = sio.load_mesh(mesh_file)
					xx = sps.quadric_curv_mean(Ks[i][0], Ks[i][1])(mesh.vertices[:, 0], mesh.vertices[:, 1])
					xx_m = np.zeros((nstep ** 2, 1))
					texture_file = cc + '_curv_mean_' + curv_type + '.gii'
					# mesh.apply_transform(mesh.principal_inertia_transform)
					os.chdir(
						"/hpc/meca/users/souani.a/curvature/Compute_curvatures/Curvatures/Noise/" + cc)
					tex = sio.load_texture(texture_file)
					err.append(mean_squared_error(xx, tex.darray))
				
				plt.plot( Sigma,err, color=couleur[s] , marker= 'o', label=str(curv_types[s]))
				print("Done for "+ cc[:-5]+ ' for '+ curv_type)
				# plt.text(Sigma[0], err[0], curv_types[s], fontsize=9)
				s += 1
			# plt.show()
			os.chdir(r"/hpc/meca/users/souani.a/curvature/Compute_curvatures/Curvatures/Figures/Mean_squarred_error_function_of_noise")
			plt.savefig("MSE function of noise for " + cc)
