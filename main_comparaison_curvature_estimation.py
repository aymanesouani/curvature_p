import matplotlib

matplotlib.use('TkAgg')

import matplotlib.pyplot as plt
import nibabel as nb
import numpy as np
import os
import slam.generate_parametric_surfaces as gp
import surf_plot as sp

if __name__ == '__main__':
	main_folder = "/hpc/meca/users/souani.a/curvature/Compute_curvatures/Quadrics"
	
	curv_types = ['Maillot', 'Rusinkiewicz', 'PetitJean', 'Taubin', 'Dong']  # , 'Peyre', 'fem', 'boix', 'bary']
	curv_types_copy = curv_types.copy()
	curv_types_copy.append('analytic')
	meshes = list()
	curv_mean = list()
	Ks = [[1, 1], [1, 0], [-1, 1]]
	
	for K in Ks:
		meshes.append('hex_quad_k1_' + str(K[0]) + '_k2_' + str(K[1]) + '_50')
		tmp_m = nb.load(os.path.join(main_folder, 'hex_quad_k1_' + str(K[0]) + '_k2_' + str(K[1]) + '_50.gii'))
		mesh_coords = tmp_m.darrays[0].data
		curv_mean.append(gp.quadric_curv_mean(K[0], K[1])(mesh_coords[:, 0], mesh_coords[:, 1]))
		print("It's done for", gp.quadric_curv_mean(K[0], K[1])(mesh_coords[:, 0], mesh_coords[:, 1]))
	
	for mesh_type, curv_ref in zip(meshes, curv_mean):
		print(mesh_type)
		surf_mesh = os.path.join(main_folder, mesh_type + '.gii')
		coords, faces = nb.load(surf_mesh).get_arrays_from_intent(nb.nifti1.intent_codes['NIFTI_INTENT_POINTSET'])[
			                0].data, \
		                nb.load(surf_mesh).get_arrays_from_intent(nb.nifti1.intent_codes['NIFTI_INTENT_TRIANGLE'])[
			                0].data
		
		vmin = np.min(curv_ref)
		vmax = np.max(curv_ref)
		fig = plt.figure(figsize=(15, 12))
		ax = fig.add_subplot(2, 3, 1)
		ax.set_title('analytic')
		sp.plot_surf_stat_map(coords, faces, stat_map=curv_ref, elev=-130, azim=30, cmap='jet', vmin=vmin,
		                            vmax=vmax,
		                            symmetric_cbar=False, ax=ax)
