import matplotlib
import nibabel
import trimesh
import slam.plot as splt
from scipy.stats import pearsonr
from sklearn.manifold import MDS

matplotlib.use('TkAgg')
import matplotlib.pyplot as plt
import nibabel as nb
import numpy as np
import os
import generate_quadric_surfaces as sps
import slam.io as sio
from matplotlib import cm
import numpy as np
from mpl_toolkits.mplot3d.axes3d import get_test_data
from mpl_toolkits.mplot3d import Axes3D
import compute_all_curvatures

if __name__ == '__main__':
	main_folder = "/hpc/meca/users/souani.a/curvature/Compute_curvatures/Noise"
	curv_types = ['Peyre', 'Rusinkiewicz', 'PetitJean', 'Taubin', 'Dong']  # , 'Peyre', 'fem', 'boix', 'bary']
	curv_types_leg = curv_types.copy()
	curv_types_leg.append('analytic')
	meshes = list()
	curv_mean = list()
	
	Ks = [[0, 1], [1, 1], [-1, 1]]
	Nstep = [10, 50, 100]  # �r exemple
	Sigma = [0.2, 0.4, 0.5, 0.55, 0.6, 0.7, 0.8]
	list_exampels = list()
	os.chdir(r"/hpc/meca/users/souani.a/curvature/Compute_curvatures/Noise")
	for i in range(len(Ks)):
		for nstep in Nstep:
			for sigma in Sigma:
				
				ks = [Ks[i]]
				X, Y, faces, Zs = sps.generate_quadric(Ks, nstep=nstep, ratio=sigma)
				Z = Zs[0]
				
				coords = np.array([X, Y, Z]).transpose()
				
				mesh = trimesh.Trimesh(faces=faces, vertices=coords, process=False)
				# mesh.show()
				
				# Estimate curvatures :
				cc = 'hex_quad_k1_' + str(Ks[i][0]) + '_k2_' + str(Ks[i][1]) + '_' + str(nstep) + '_noise_' + str(
					int(sigma * 1000))
				# os.chdir(r"/hpc/meca/users/souani.a/curvature/Compute_curvatures/Noise")
				# os.mkdir(cc )
				os.chdir(r"/hpc/meca/users/souani.a/curvature/Compute_curvatures/Noise/"+cc)
				sio.write_mesh(mesh, cc+'.gii')
				list_exampels.append(cc)
				print("done for " + cc)
				data = list()
				for curv_type in curv_types:
					mesh_file = 'hex_quad_k1_' + str(Ks[i][0]) + '_k2_' + str(Ks[i][1]) + '_' + str(nstep) + '.gii'
					texture_file = cc + '_curv_mean_' + curv_type + '.gii'
					os.chdir("/hpc/meca/users/souani.a/curvature/Compute_curvatures/Quadrics")
					mesh = sio.load_mesh(mesh_file)
					# mesh.apply_transform(mesh.principal_inertia_transform)
					
					os.chdir(
						"/hpc/meca/users/souani.a/curvature/Compute_curvatures/Noise/" + cc)
					tex = sio.load_texture(texture_file)
					
					# splt.pyglet_plot(mesh, tex.darray, plot_colormap=True)
					# splt.pyglet_plot(mesh, tex.darray, 'hot', plot_colormap=True)
					# print(tex.darray)
					data.append(tex.darray)
				xx = sps.quadric_curv_mean(Ks[i][0], Ks[i][1])(mesh.vertices[:, 0], mesh.vertices[:, 1])
				xx = np.reshape(xx, (1, nstep ** 2))
				data = np.array(data)
				data_m = np.zeros((6, nstep ** 2))
				data_m[:5, :] = data[:, :, 0]
				data_m[5, :] = xx
				
				
				data_m = np.array(data_m)
				embedding = MDS(n_components=2)
				data_t = embedding.fit_transform(data_m[:, :])
				# print(data_t)
				fig = plt.figure()
				ax = fig.add_subplot(111)
				fig.suptitle(cc)
				plt.xlabel('First component')
				plt.ylabel('Second component')
				x, y = data_t[:, 0], data_t[:, 1]
				plt.plot(x, y, 'ro')
				for k in range(6):
					plt.scatter(x[k], y[k], marker='x', color='red')
					plt.text(x[k], y[k], curv_types_leg[k], fontsize=9)
				plt.savefig(cc + "+analytique")
	
	# sio.write_mesh(mesh, 'hex_quad_k1_'+str(Ks[0][0])+'_k2_'+str(Ks[0][1])+'_'+str(nstep)+'.gii')
	# print(mesh.vertices == coords)
	# os.chdir(r"/hpc/meca/users/souani.a/curvature/Compute_curvatures/Noise")
	"""for i in range(len(list_exampels)):
		os.chdir(r"/hpc/meca/users/souani.a/curvature/Compute_curvatures/Noise/"+list_exampels[i])
		file_out = list_exampels[i]
		output_folder = '/hpc/meca/users/souani.a/curvature/Compute_curvatures/Noise/'+list_exampels[i]
		compute_all_curvatures.matlab_curvatures(list_exampels[i] + '.gii', file_out, output_folder)
		print("Done for "+ list_exampels[i])
	print("It comes to an end")"""

	for nstep in Nstep:
		for i in range(len(Ks)):
			for sigma in Sigma:
				meshes.append('hex_quad_k1_'+str(Ks[i][0])+'_k2_'+str(Ks[i][1])+'_'+str(nstep)+'_noise_'+str(int(sigma * 1000)))
				os.chdir(r"/hpc/meca/users/souani.a/curvature/Compute_curvatures/Quadrics")
				tmp_m = sio.load_mesh('hex_quad_k1_'+str(Ks[i][0])+'_k2_'+str(Ks[i][1])+'_'+str(nstep)+'.gii')
				mesh_coords = tmp_m.vertices
				curv_mean.append(sps.quadric_curv_mean(Ks[i][0], Ks[i][1])(mesh_coords[:, 0], mesh_coords[:, 1]))
				
	