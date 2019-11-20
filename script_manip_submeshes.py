import os

from sklearn.manifold import MDS
import matplotlib.pyplot as plt
import Rusinkei.CalculCurvature as CC
import slam.io as sio
import trimesh
import real_mesh_compute_all_curvatures
import compute_all_curvatures
import slam.plot as splt
import numpy as np
import generate_quadric_surfaces as gqs




if __name__ == "__main__" :
	ccc = ["sub_mesh0", "sub_mesh1", "sub_mesh2", "sub_mesh3", "sub_mesh4"]
	curv_types = ["Dong", "PetitJean", "Peyre", "Rusinkiewicz", "Taubin"]
	mm = "lhwhite"
	i = 0
	for cc in ccc:
		
		os.chdir(r"/hpc/meca/users/souani.a/curvature")
		data = list()
		mesh = sio.load_mesh(cc+".gii")
		for curv_type in curv_types:
			print(curv_type)
			os.chdir(r"/hpc/meca/users/souani.a/data/sub_m/sub_m"+str(i))
			texture_mean = sio.load_texture(cc+"_curv_mean_"+ curv_type+".gii")
			print(np.shape(texture_mean.darray))
			data.append(texture_mean.darray)
			# splt.pyglet_plot(mesh, texture_mean.darray, color_map='jet')
			texture_gauss = sio.load_texture(cc + "_curv_gauss_" + curv_type + ".gii")
			# splt.pyglet_plot(mesh, texture_gauss.darray, color_map='jet')
		#data = np.array(data)
		#print(data)
		data = np.array(data)
		print(np.shape(data))
		embedding = MDS(n_components=2)
		print(np.shape(data[:, :, 0]))
		data_t = embedding.fit_transform(data[:, :, 0])
		print(np.shape(data_t))
		# print(data_t)
		fig = plt.figure()
		ax = fig.add_subplot(111)
		fig.suptitle(cc)
		plt.xlabel('First component')
		plt.ylabel('Second component')
		print(np.shape(data_t))
		x, y = data_t[:, 0], data_t[:, 1]
		plt.plot(x, y, 'ro')
		for k in range(5):
			plt.text(x[k], y[k], curv_types[k], fontsize=9)
		print("surface = ", mesh.area_faces.sum())
		os.chdir("/hpc/meca/users/souani.a/curvature/sub_meshes/Figures")
		plt.savefig(cc + " isomaps de surface "+ str(int(mesh.area_faces.sum())))
		plt.show()
		i+= 1
	print("Done for " + mm)
