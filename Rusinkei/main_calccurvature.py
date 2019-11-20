import trimesh
from trimesh.curvature import discrete_gaussian_curvature_measure, discrete_mean_curvature_measure, sphere_ball_intersection
import matplotlib.pyplot as plt
import numpy as np
import Rusinkei.CalculCurvature as CC
import generate_quadric_surfaces as gp
import matplotlib as mpl

if __name__ == '__main__':
	
	# Generate a quadric
	X, Y, faces, Zs = gp.generate_quadric([[-1, 1]], nstep=100)
	Z=Zs[0]
	coords = np.array([X, Y, Z]).transpose()
	mesh = trimesh.Trimesh(faces=faces, vertices=coords, process=False)
	
	
	# Show th sphere
	mesh.show()
	
	# Calculate Rusinkiewicz estimation of mean and gauss curvatures
	PrincipalCurvatures, PrincipalDir1, PrincipalDir2 = CC.GetCurvaturesAndDerivatives(mesh)
	gaussian_curv = PrincipalCurvatures[0, :] * PrincipalCurvatures[1, :]
	mean_curv = 0.5 * (PrincipalCurvatures[0, :] + PrincipalCurvatures[1, :])
	
	# Plot mean curvature
	vect_col_map = \
		trimesh.visual.color.interpolate(mean_curv, color_map='jet')
	
	if mean_curv.shape[0] == mesh.vertices.shape[0]:
		mesh.visual.vertex_colors = vect_col_map
	elif mean_curv.shape[0] == mesh.faces.shape[0]:
		mesh.visual.face_colors = vect_col_map
	mesh.show( background=[0, 0, 0, 255])
	
	
	mpl.rcParams.update({'font.size': 20})
	# fig = plt.figure(figsize=(8, 2))
	# ax1 = fig.add_axes([0.05, 0.80, 0.9, 0.15])
	fig, ax = plt.subplots(1, 1)
	ax.set_title("Colormap")
	
	
	vmin = np.min(mean_curv)
	vmax = np.max(mean_curv)
	norm = mpl.colors.Normalize(vmin=vmin, vmax=vmax)
	
	mpl.colorbar.ColorbarBase(ax, norm=norm, orientation = 'horizontal')
	fig.set_size_inches(18, 4)
	plt.show()
	
	# PLot Gauss curvature
	vect_col_map = \
		trimesh.visual.color.interpolate(gaussian_curv, color_map='jet')
	if gaussian_curv.shape[0] == mesh.vertices.shape[0]:
		mesh.visual.vertex_colors = vect_col_map
	elif gaussian_curv.shape[0] == mesh.faces.shape[0]:
		mesh.visual.face_colors = vect_col_map
	mesh.show(background=[0, 0, 0, 255])

	mpl.rcParams.update({'font.size': 20})
	# fig = plt.figure(figsize=(8, 2))
	# ax1 = fig.add_axes([0.05, 0.80, 0.9, 0.15])
	fig, ax = plt.subplots(1, 1)
	ax.set_title("Colormap")
	
	vmin = np.min(gaussian_curv)
	vmax = np.max(gaussian_curv)
	norm = mpl.colors.Normalize(vmin=vmin, vmax=vmax)
	
	mpl.colorbar.ColorbarBase(ax, norm=norm, orientation='horizontal')
	fig.set_size_inches(18, 4)
	plt.show()
	# Comparaison between estimated and  real curvatures
	
	
"""	print("mean error ", np.linalg.norm(mean - mean_curv))
	print("gauss error ", np.linalg.norm(gauss - gaussian_curv))"""
