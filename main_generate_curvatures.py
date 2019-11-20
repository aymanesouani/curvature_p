import random

from sklearn.manifold import MDS
from trimesh import Trimesh

import compute_all_curvatures
import generate_quadric_surfaces as gqs
import os
import slam.io as sio
import numpy as np
import matplotlib.pyplot as plt
import trimesh


def f(x):
	return 0.5 * x * np.sqrt(1 + 4 * x ** 2) + 0.25 * np.log(2 * x + np.sqrt(1 + 4 * x ** 2))


if __name__ == '__main__':
	os.chdir(r"/hpc/meca/users/souani.a/curvature/Compute_curvatures/Curvatures")
	list_nbr_steps = [46, 32, 22, 19, 15]
	list_surfaces = [1456.5160195140602, 642.8849596346015, 331.0259595916111, 246.50582433706595, 140.8828640008205]
	list_densities = np.array(list_nbr_steps) ** 2 / np.array(list_surfaces)
	print(list_densities)
	Ks = [[0, 1], [1, 1], [-1, 1]]
	list_ax_ay = list()
	liste_couples = np.arange(1, 10, 0.05)
	# print (liste_couples)
	eps = 1e-1
	dico = {1: [(8.85, 6.25), (9.9, 3.8), (9.9, 2.6), (9.65, 2.2), (9.75, 1.5)], 2 : [(9.65, 3.4), (9.75, 1.55), (8.65, 1.), (7.4, 1.), (5.45, 1.)], 3: [(9.65, 3.4), (9.6, 1.6), (8.45, 1.05), (7.4, 1.0), (5.45, 1.)]}
	minimum = (liste_couples[0], liste_couples[0])
	curv_types = ['Peyre', 'Rusinkiewicz', 'PetitJean', 'Taubin', 'Dong']
	curv_types_leg = curv_types.copy()
	curv_types_leg.append("Analytique")
	Sigma = [0.2, 0.4, 0.5, 0.55, 0.6, 0.7, 0.8]
	
	for ks in Ks:
		for i in range(len(list_nbr_steps)):
			for sigma in Sigma:
			
					# Generate the quadric
					# X, Y, faces, Zs = gqs.generate_quadric([ks], nstep=list_nbr_steps[i], ax=list_ax_ay[i][0ay=list_ax_ay[i][1])
					# Z = Zs[0]
					# coords = np.array([X, Y, Z]).transpose()
					# mesh = trimesh.Trimesh(faces=faces, vertices=coords, process=False)
					# mesh.show()
					# print("mesh surface = ", mesh.area, "\n surface needed = ", list_surfaces[i])
					
					os.chdir(r"/hpc/meca/users/souani.a/curvature/Compute_curvatures/Curvatures/0Noise")
					cc = 'hex_quad_k1_' + str(ks[0]) + '_k2_' + str(ks[1]) + '_' + str(
						list_nbr_steps[i])
					mesh = sio.load_mesh(cc+'.gii')
					data = list()
					xx = gqs.quadric_curv_mean(ks[0], ks[1])(mesh.vertices[:, 0], mesh.vertices[:, 1])
					print("-------------------------------" + cc + '_noise_' + str(int(sigma * 1000)))
					for  curv_type in curv_types :
						print(curv_type)
						texture_file = cc + '_noise_' + str(int(sigma * 1000))+ '_curv_mean_' + curv_type + '.gii'
						
						# mesh.apply_transform(mesh.principal_inertia_transform)
						
						os.chdir(
							"/hpc/meca/users/souani.a/curvature/Compute_curvatures/Curvatures/Noise/" + cc + '_noise_' + str(int(sigma * 1000)))
						tex = sio.load_texture(texture_file)
						data.append(tex.darray)
					data = np.array(data)
					xx = np.reshape(xx, (1, list_nbr_steps[i] ** 2))
					data_m = np.zeros((6, list_nbr_steps[i] ** 2))
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
					plt.savefig(cc + " + analytique")
					
					# j+=1
		
	
	'''Sigma = [0.2, 0.4, 0.5, 0.55, 0.6, 0.7, 0.8]
	for sigma in Sigma:
		k = 1
		for ks in Ks:
			list_ax_ay = dico[k]
			for i in range(len(list_nbr_steps)):
				os.chdir(r"/hpc/meca/users/souani.a/curvature/Compute_curvatures/Curvatures/Noise")
				X, Y, faces, Zs = gqs.generate_quadric([ks], nstep=list_nbr_steps[i], ax=list_ax_ay[i][0],
				                                       ay=list_ax_ay[i][1], ratio= sigma)
				Z = Zs[0]
				coords = np.array([X, Y, Z]).transpose()
				mesh = trimesh.Trimesh(faces=faces, vertices=coords, process=False)
				# mesh.show()
				cc = 'hex_quad_k1_'+str(ks[0])+'_k2_'+ str(ks[1])+ '_'+str(list_nbr_steps[i])+ '_noise_' + str(
					int(sigma * 1000))
				# sio.write_mesh(mesh, cc + '.gii')
				print("mesh surface = ", mesh.area, "\n surface needed = ", list_surfaces[i], "\n difference =", mesh.area - list_surfaces[i])
				print("density =", mesh.density)
				
				# Compute estimations via Matlab
				# os.mkdir(cc)
				# file_out = cc
				# output_folder = '/hpc/meca/users/souani.a/curvature/Compute_curvatures/Curvatures/Noise/' + cc
				# compute_all_curvatures.matlab_curvatures(cc+ '.gii', file_out, output_folder)
				print("Done for " + cc)
			
			k += 1
	
	print("It comes to an end")'''
