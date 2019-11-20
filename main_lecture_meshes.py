import os
import nibabel
import numpy as np
from matplotlib.collections import LineCollection
from scipy.stats import pearsonr
from sklearn.decomposition import PCA
from sklearn.manifold import MDS
from sklearn.metrics import euclidean_distances
from trimesh.visual import color
from time import time
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
from matplotlib.ticker import NullFormatter
from sklearn import manifold



def corr_vect(a, b) :
	a = np.asarray(a)
	b = np.asarray(b)
	a_mean = a.mean()
	b_mean = b.mean()
	a_centered = a - a_mean
	b_centered = b - b_mean
	r = np.dot(a_centered / np.linalg.norm(a_centered), b_centered / np.linalg.norm(b_centered))
	r = max(min(r, 1.0), -1.0)
	return r


def plot_embedding(X, title=None):
	x_min, x_max = np.min(X, 0), np.max(X, 0)
	X = (X - x_min) / (x_max - x_min)
	
	plt.figure()
	ax = plt.subplot(111)
	for i in range(X.shape[0]):
		plt.text(X[i, 0], X[i, 1], str(y[i]),
		         color=plt.cm.Set1(y[i] / 10.),
		         fontdict={'weight': 'bold', 'size': 9})
	
	if hasattr(offsetbox, 'AnnotationBbox'):
		# only print thumbnails with matplotlib > 1.0
		shown_images = np.array([[1., 1.]])  # just something big
		for i in range(X.shape[0]):
			dist = np.sum((X[i] - shown_images) ** 2, 1)
			if np.min(dist) < 4e-3:
				# don't show points that are too close
				continue
			shown_images = np.r_[shown_images, [X[i]]]
			imagebox = offsetbox.AnnotationBbox(
				offsetbox.OffsetImage(digits.images[i], cmap=plt.cm.gray_r),
				X[i])
			ax.add_artist(imagebox)
	plt.xticks([]), plt.yticks([])
	if title is not None:
		plt.title(title)

if __name__ == '__main__':
	os.chdir(r"/hpc/meca/users/bohi.a/Codes/PycharmProjects/CurvatureEstimation")
	os.chdir(r"/hpc/meca/users/souani.a/curvature/Compute_curvatures")
	# 'hex_quad_k1_0_k2_0_10', 'hex_quad_k1_0_k2_0_50',
	list_exampels = [ 'hex_quad_k1_0_k2_1_10', 'hex_quad_k1_0_k2_1_50',
	                 'hex_quad_k1_0_k2_-1_10', 'hex_quad_k1_0_k2_-1_50', 'hex_quad_k1_1_k2_0_10',
	                 'hex_quad_k1_1_k2_0_50', 'hex_quad_k1_1_k2_1_10', 'hex_quad_k1_1_k2_1_50', 'hex_quad_k1_1_k2_1_50',
	                 'hex_quad_k1_1_k2_-1_50', 'hex_quad_k1_-1_k2_0_10', 'hex_quad_k1_-1_k2_0_50',
	                 'hex_quad_k1_-1_k2_1_10', 'hex_quad_k1_-1_k2_1_50', 'hex_quad_k1_-1_k2_-1_10',
	                 'hex_quad_k1_-1_k2_-1_50']
	print(len(list_exampels))
	# for x in list_exampels :
	# print("working on quadric " + x)
	# path_file_mesh =  x+'.gii'
	# output_file = '/hpc/meca/users/souani.a/curvature/Compute_curvatures'
	# file_out = x
	# compute_all_curvatures.matlab_curvatures(path_file_mesh, file_out, output_file)
	# print("comparaison process finished")
	curv_types = ['Peyre', 'Rusinkiewicz', 'PetitJean', 'Taubin', 'Dong']  # Maillot / Peyre ??
	for quadric in list_exampels:
		data = list()
		mx_distances_quadric = np.zeros((5, 5))
		mx_corr_quadric = np.zeros((5, 5))
		x_corr_quadric = np.zeros((5, 5))
		
		for i in range (len(curv_types)):
			texture_moy_1 = nibabel.load(os.path.join("/hpc/meca/users/souani.a/curvature/Compute_curvatures",
			                                          quadric + '_curv_mean_' + curv_types[i] + '.gii'))
			curv_vals_1 = texture_moy_1.darrays[0].data.squeeze()
			if curv_types[i] is 'PetitJean':
				curv_vals_1 = curv_vals_1 * 2
			data.append(curv_vals_1)
			for j in range (len(curv_types)):
			
				"texture_gauss = sio.load_texture(quadric+'_curv_gauss_'+curv_type+'.gii')"
				
				texture_moy_2 = nibabel.load(os.path.join("/hpc/meca/users/souani.a/curvature/Compute_curvatures", quadric + '_curv_mean_' + curv_types[j] + '.gii'))
				
				curv_vals_2 = texture_moy_2.darrays[0].data.squeeze()
				if curv_types[j] is 'PetitJean':
					curv_vals_2 = 2 * curv_vals_2
				mx_distances_quadric[i, j] = np.linalg.norm(curv_vals_1 - curv_vals_2)
				mx_corr_quadric[i, j] = corr_vect(curv_vals_1, curv_vals_2)
				x_corr_quadric[i, j] = pearsonr(curv_vals_1, curv_vals_2)[0]
		print(quadric +"\n" , mx_distances_quadric)
		print(quadric + "\n", mx_corr_quadric)
		print(quadric + "\n", x_corr_quadric)
		" Generate  data "
		data = np.array(data)
		embedding = MDS(n_components=2)
		data_transformed = embedding.fit_transform(data)
		print(data_transformed)
		print("shape data transformed", np.shape(data_transformed))
		print("Matrice de distance MDS", embedding.metric)
		
		fig = plt.figure()
		ax = fig.add_subplot(111)
		fig.suptitle(quadric)
		plt.xlabel('First component')
		plt.ylabel('Second component')
		x, y = data_transformed[:, 0], data_transformed[:, 1]
		plt.plot(x, y, 'ro')
		for i in range(5):
			plt.scatter(x[i], y[i], marker='x', color='red')
			plt.text(x[i], y[i], curv_types[i], fontsize=9)
		plt.savefig(quadric)
		
		# Calcul with PCA
		pca = PCA(n_components=3)
		coords = pca.fit_transform(data)
		# pca.fit(data)
		
		eigen_values = pca.explained_variance_
		print("eigen values ", eigen_values)
		print("First component inertie",  100 * eigen_values[0] / sum(eigen_values))
		print("Second component inertie", 100 * eigen_values[1] / sum(eigen_values))
		print("Third component inertie",  100 * eigen_values[2] / sum(eigen_values))
		
		# ax.plot(x, y, 'ro')
		
		# ax.set_xlabel('First component')
		# ax.set_ylabel('Second component')
		# ax.set_zlabel('Third component')
		# plt.title(quadric+" methods projection")
		"""for i in range(5):  # plot each point + it's index as text above
			 # fig.scatter(x[i], y[i],, color='b')
			 fig.text(x[i], y[i], '%s' % curv_types[i], size=20, zorder=1,
			        color='k')
		
		# Add noise to the similarities
		# don't add noise to the similarities
		
		
		"""
		
		# plt.savefig(quadric+'2d')
		"""
for quadric in list_exampels:
		data = list()
		for i in range(len(curv_types)):
			"texture_gauss = sio.load_texture(quadric+'_curv_gauss_'+curv_type+'.gii')"
			texture_moy_1 = nibabel.load(os.path.join("/hpc/meca/users/souani.a/curvature/Compute_curvatures",
			quadric + '_curv_mean_' + curv_types[i] + '.gii'))
			curv_vals_1 = texture_moy_1.darrays[0].data.squeeze()
			data.append(curv_vals_1)
		data = np.array(data)
		print(np.shape(data))
		X = data
		y = np.transpose(np.array([1, 2, 3, 4, 5]))
		n_features, n_samples = np.shape(X)
		n_neighbors = 2
		fig = plt.figure(figsize=(15, 8))
		plt.suptitle("Manifold Learning with %i points, %i neighbors"
		             % (1000, n_neighbors), fontsize=14)
		
		# ----------------------------------------------------------------------
		# Scale and visualize the embedding vectors
		
		
		
		# ----------------------------------------------------------------------
		# Plot images of the digits
		t0 = time()
		mds = manifold.MDS(n_components = 2, max_iter=100, n_init=1)
		Y = mds.fit_transform(X)
		t1 = time()
		print("MDS: %.2g sec" % (t1 - t0))
		ax = fig.add_subplot(258)
		plt.scatter(Y[:, 0], Y[:, 1], c=color, cmap=plt.cm.Spectral)
		plt.title("MDS (%.2g sec)" % (t1 - t0))
		ax.xaxis.set_major_formatter(NullFormatter())
		ax.yaxis.set_major_formatter(NullFormatter())
		plt.axis('tight')
		plt.show()"""