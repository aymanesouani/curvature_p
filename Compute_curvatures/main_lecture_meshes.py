
import trimesh
import slam.io as sio
import slam.plot as splt
import nibabel as nb
import os
import numpy as np

if __name__ == '__main__':
	fd = '/hpc/meca/users/bohi.a/Data/Models/week23_ref_surf/B0.txt'
	coords = np.loadtxt(fd, skiprows=1, max_rows=50943)
	faces = np.loadtxt(fd, skiprows=50945, dtype=np.int)
	mesh = trimesh.Trimesh(faces=faces - 1, vertices=coords, process=False)
	mesh.show()
	os.chdir(r"/hpc/meca/users/souani.a/data/OASIS")
	mesh1 = sio.load_mesh('OAS1_0277_Rwhite.gii')
	mesh1.show()
	# curv_gifti = nb.gifti.read('hex_quad_k1_0_k2_10_curv_mean_Dong.gii')
	# curv_tex = curv_gifti.darrays[0].data.squeeze()
	
	# splt.pyglet_plot(mesh1, curv_tex)
	"""mesh1 = sio.load_mesh('example_1_curv_mean_Dong.gii')
	mesh1.show()
	mesh1 = sio.load_mesh('example_1_curv_gauss_PetitJean.gii')
	mesh1.show()
	mesh1 = sio.load_mesh('example_1_curv_mean_PetitJean.gii')
	mesh1.show()
	mesh1 = sio.load_mesh('example_1_curv_gauss_Peyre.gii')
	mesh1.show()
	mesh1 = sio.load_mesh('example_1_curv_mean_Peyre.gii')
	mesh1.show()
	mesh1 = sio.load_mesh('example_1_curv_gauss_Taubin.gii')
	mesh1.show()
	mesh1 = sio.load_mesh('example_1_curv_mean_Taubin.gii')
	mesh1.show()
	mesh1 = sio.load_mesh('example_1_curv_gauss_Rusinkiewicz.gii')
	mesh1.show()
	mesh1 = sio.load_mesh('example_1_curv_mean_Rusinkiewicz.gii')
	mesh1.show()"""

