import os

import compute_all_curvatures

if __name__ == '__main__':
	print("Début script")
	list_exampels = [ 'hex_quad_k1_0_k2_-1_100',
	                 'hex_quad_k1_1_k2_-1_100',
	                 'hex_quad_k1_1_k2_0_100',
	                 'hex_quad_k1_1_k2_1_100',
	                 'hex_quad_k1_-1_k2_0_100',
	                 'hex_quad_k1_-1_k2_1_100',
	                 'hex_quad_k1_-1_k2_-1_100']
	os.chdir(r"/hpc/meca/users/souani.a/curvature/Compute_curvatures/Quadrics")
	
	for i in range (len(list_exampels)):
		os.chdir(r"/hpc/meca/users/souani.a/curvature/Compute_curvatures/Quadrics")
		file_out = list_exampels[i]
		output_file = '/hpc/meca/users/souani.a/curvature/Compute_curvatures'
		compute_all_curvatures.matlab_curvatures(list_exampels[i] + '.gii', file_out, output_file)
		print("Done for "+ list_exampels[i])
	print("It comes to an end")
