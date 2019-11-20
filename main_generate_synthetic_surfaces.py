import numpy as np
import slam.io as sio
import os

if __name__ == '__main__':
	list_str_meshes = ['sub_mesh0', 'sub_mesh1', 'sub_mesh2', 'sub_mesh3', 'sub_mesh4']
	data_path = "/hpc/meca/users/souani.a/curvature/"
	list_surfaces, list_densities, list_nbr_vertex = list(), list(), list()
	os.chdir(data_path)
	for str_mesh in list_str_meshes:
		mesh = sio.load_mesh(str_mesh + '.gii')
		# mesh.show()
		print(str_mesh, ' surface : ', mesh.area_faces.sum())
		list_surfaces.append(mesh.area_faces.sum())
		print(str_mesh, ' vertices nbr : ', np.shape(mesh.vertices)[0])
		list_nbr_vertex.append(np.shape(mesh.vertices)[0])
		print(str_mesh, ' density : ', np.shape(mesh.vertices)[0] / mesh.area_faces.sum())
		list_densities.append(np.shape(mesh.vertices)[0] / mesh.area_faces.sum())
	# Comparaison d'histogrammes
	
