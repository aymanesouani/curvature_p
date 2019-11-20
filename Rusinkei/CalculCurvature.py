import time
import matplotlib.pyplot as plt
import numpy as np
import trimesh
import matplotlib as mpl

#import slam.io as sio
#from slam.plot import pyglet_plot


def ProjectCurvatureTensor(uf, vf, nf, old_ku, old_kuv, old_kv, up, vp):
	r_new_u, r_new_v = RotateCoordinateSystem(up, vp, nf)
	OldTensor = np.array([[old_ku, old_kuv], [old_kuv, old_kv]])
	u1 = np.dot(r_new_u, uf)
	v1 = np.dot(r_new_u, vf)
	u2 = np.dot(r_new_v, uf)
	v2 = np.dot(r_new_v, vf)
	
	new_ku = np.dot(np.array([u1, v1]), np.dot(OldTensor, np.transpose(np.array([u1, v1]))))
	new_kuv = np.dot(np.array([u1, v1]), np.dot(OldTensor, np.transpose(np.array([u2, v2]))))
	new_kv = np.dot(np.array([u2, v2]), np.dot(OldTensor, np.transpose(np.array([u2, v2]))))
	return new_ku, new_kuv, new_kv


def CalcCurvature(FV, VertexNormals, FaceNormals, Avertex, Acorner, up, vp):
	"""CalcFaceCurvature recives a list of vertices and faces in FV structure
	and the normal at each vertex and calculates the second fundemental
	matrix and the curvature using least squares
	
	INPUT :
	FV - face-vertex data structure containing a list of vertices and a list of faces
	VertexNoRMALS - n*3 matrix ( n = number of vertices ) containing the normal at each vertex
	FaceNormals - m*3 matrix ( m = number of faces ) containing the normal of each face
	
	OUTPOUT
	FaceSFM - an m*1 cell matrix second fundemental
	VertexSFM - an n*w cell matrix second fundementel
	wfp - corner voronoi weights """
	
	print(" Calculating curvature tensors ..... Please wait")
	"Matrix of each face at each cell"
	FaceSFM, VertexSFM = list(), list()
	for i in range(FV.faces.shape[0]):
		FaceSFM.append([[0, 0], [0, 0]])
	for i in range(FV.vertices.shape[0]):
		VertexSFM.append([[0, 0], [0, 0]])
	Kn = np.zeros((1, FV.faces.shape[0]))
	
	" Get all the edge vectors "
	e0 = FV.vertices[FV.faces[:, 2], :] - FV.vertices[FV.faces[:, 1], :]
	e1 = FV.vertices[FV.faces[:, 0], :] - FV.vertices[FV.faces[:, 2], :]
	e2 = FV.vertices[FV.faces[:, 1], :] - FV.vertices[FV.faces[:, 0], :]
	
	" Normalize edge vectors "
	e0_norm = normr(e0)
	# e1_norm = normr(e1)
	# e2_norm = normr(e2)
	
	wfp = np.array(np.zeros((FV.faces.shape[0], 3)))
	for i in range(FV.faces.shape[0]):
		"Calculate Curvature Per Face"
		"set face coordinate frame"
		nf = FaceNormals[i, :]
		t = e0_norm[i, :]
		B = np.cross(nf, t)
		B = B / (np.linalg.norm(B))
		
		"extract relevant normals in face vertices"
		n0 = VertexNormals[FV.faces[i][0], :]
		n1 = VertexNormals[FV.faces[i][1], :]
		n2 = VertexNormals[FV.faces[i][2], :]
		
		" solve least squares problem of th form Ax=b "
		A = np.array([[np.dot(e0[i, :], t), np.dot(e0[i, :], B), 0], [0, np.dot(e0[i, :], t), np.dot(e0[i, :], B)],
		              [np.dot(e1[i, :], t), np.dot(e1[i, :], B), 0], [0, np.dot(e1[i, :], t), np.dot(e1[i, :], B)],
		              [np.dot(e2[i, :], t), np.dot(e2[i, :], B), 0], [0, np.dot(e2[i, :], t), np.dot(e2[i, :], B)]])
		
		b = np.array(
			[np.dot(n2 - n1, t), np.dot(n2 - n1, B), np.dot(n0 - n2, t), np.dot(n0 - n2, B), np.dot(n1 - n0, t),
			 np.dot(n1 - n0, B)])
		
		"Resolving by least mean square method because A is not a square matrix"
		
		x = np.linalg.lstsq(A, b, None)
		
		FaceSFM[i] = np.array([[x[0][0], x[0][1]], [x[0][1], x[0][2]]])
		Kn[0][i] = np.dot(np.array([1, 0]), np.dot(FaceSFM[i], np.array([[1.], [0.]])))
		"""
		Calculate curvature per vertex
		Calculate voronoi weights
		"""
		wfp[i][0] = Acorner[i][0] / Avertex[FV.faces[i][0]]
		wfp[i][1] = Acorner[i][1] / Avertex[FV.faces[i][1]]
		wfp[i][2] = Acorner[i][2] / Avertex[FV.faces[i][2]]
		
		"Calculate new coordinate system and project the tensor"
		
		for j in range(3):
			new_ku, new_kuv, new_kv = ProjectCurvatureTensor(t, B, nf, x[0][0], x[0][1], x[0][2],
			                                                 up[FV.faces[i][j], :], vp[FV.faces[i][j], :])
			VertexSFM[FV.faces[i][j]] += np.dot(wfp[i][j], np.array([[new_ku, new_kuv], [new_kuv, new_kv]]))
	
	print('Finished Calculating curvature tensors')
	
	return FaceSFM, VertexSFM, wfp


"""def CalculDerivative(FV, FaceNormals, PrincipalCurvatures, up, vp, wfp):
	
	"CalcCurvatureDerivative recives a list of vertices and faces in FV structure
	"and the curvature tensor at each vertex and calculate the curvature
	"derivative matrix 2*2*2 VertexCmatrix using least squares
	
	
	
	"Input:
	"FV - face-vertex data structure containing a list of vertices and a list of faces
	"VertexSFM - n*(2*2) cell array (n= number of vertices) containing the second funcdmental tensor at each vertex
	"FaceSFM - m*(2*2) cell array (m=number of faces) containing the second fundmental tensor at each face
	
	
	
	"Output
	"FaceCmatrix - an  m*(2*2*2) cell matrix (m = number of faces) second fundemental derivative tensor at each face
	"VertexCmatrix - an  n*(2*2*2) cell matrix (n = number of vertices) second fundemental derivative tensor at each vertex
	
	
	print ('Calculating C Tensors ... Please wait')
	
	FaceCMatrix = np.zeros((FV.faces.shape[0], 4))
	VertexCMatrix = np.zeros((FV.vertices.shape[0], 4))
	new_CMatrix = np.zeros((1, 4))
	
	return None
	"""


def GetCurvaturesAndDerivatives(FV):
	FaceNormals = CalcFaceNormals(FV)
	(VertexNormals, Avertex, Acorner, up, vp) = CalcVertexNormals(FV, FaceNormals)
	(FaceSFM, VertexSFM, wfp) = CalcCurvature(FV, VertexNormals, FaceNormals, Avertex, Acorner, up, vp)
	[PrincipalCurvature, PrincipalDi1, PrincipalDi2] = getPrincipalCurvatures(FV, VertexSFM, up, vp)
	return PrincipalCurvature, PrincipalDi1, PrincipalDi2


def CalcFaceNormals(FV):
	"""
	
	:param FV: face vertex data structure containing a list of vertices and a list of faces
	:return:
	"""
	
	"Get all edge vectors"
	e0 = FV.vertices[FV.faces[:, 2], :] - FV.vertices[FV.faces[:, 1], :]
	e1 = FV.vertices[FV.faces[:, 0], :] - FV.vertices[FV.faces[:, 2], :]
	
	"Calculate and return normal of face"
	"FaceNormals = np.cross(e0, e1)"
	return normr(np.cross(e0, e1))


def normr(X):
	"""
	Returns a matrix with the same size where each row normalized to a vector length of 1
	"""
	
	if len(np.shape(X)) == 1:
		return X / np.abs(X)
	else:
		a = np.shape(X)[1]
		b = np.shape(X)[0]
		return np.dot(np.reshape(np.transpose(np.sqrt(1 / somme_colonnes(np.transpose(X ** 2)))), (b, 1)),
		              np.ones((1, a))) * X


def CalcVertexNormals(FV, N):
	print("Calculating vertex normals .... Please wait")
	
	"Get all the edge vectors"
	e0 = np.array(FV.vertices[FV.faces[:, 2], :] - FV.vertices[FV.faces[:, 1], :])
	e1 = np.array(FV.vertices[FV.faces[:, 0], :] - FV.vertices[FV.faces[:, 2], :])
	e2 = np.array(FV.vertices[FV.faces[:, 1], :] - FV.vertices[FV.faces[:, 0], :])
	
	"Normalize edge vectors "
	e0_norm = normr(e0)
	e1_norm = normr(e1)
	e2_norm = normr(e2)
	
	de0 = np.sqrt((e0[:, 0]) ** 2 + (e0[:, 1]) ** 2 + (e0[:, 2]) ** 2)
	de1 = np.sqrt((e1[:, 0]) ** 2 + (e1[:, 1]) ** 2 + (e1[:, 2]) ** 2)
	de2 = np.sqrt((e2[:, 0]) ** 2 + (e2[:, 1]) ** 2 + (e2[:, 2]) ** 2)
	l2 = np.array([de0 ** 2, de1 ** 2, de2 ** 2])
	l2 = np.transpose(l2)
	"""
	using ew to calulate the cot of the angles for the voronoi area calculation
	ew is the triangle barycenter. We check later if it's inside or outside the triangle
	"""
	ew = np.array([l2[:, 0] * (l2[:, 1] + l2[:, 2] - l2[:, 0]), l2[:, 1] * (l2[:, 2] + l2[:, 0] - l2[:, 1]),
	               l2[:, 2] * (l2[:, 0] + l2[:, 1] - l2[:, 2])])
	s = (de0 + de1 + de2) / 2
	
	"Af - face area vector"
	Af = np.sqrt(s * (s - de0) * (s - de1) * (s - de2))
	
	"herons formula for triangle area, could have also used  0.5 * norm(cross(e0,e1)) "
	"Calc weights"
	Acorner = np.zeros((np.shape(FV.faces)[0], 3))
	Avertex = np.zeros((np.shape(FV.vertices)[0], 1))
	
	"Calcul vertices normals"
	VertexNormals, up, vp = np.zeros((np.shape(FV.vertices)[0], 3)), np.zeros((np.shape(FV.vertices)[0], 3)), np.zeros(
		(np.shape(FV.vertices)[0], 3))
	
	for i in range(np.shape(FV.faces)[0]):
		wfv1 = Af[i] / ((de1[i] ** 2) * (de2[i] ** 2))
		wfv2 = Af[i] / ((de0[i] ** 2) * (de2[i] ** 2))
		wfv3 = Af[i] / ((de1[i] ** 2) * (de0[i] ** 2))
		
		VertexNormals[FV.faces[i][0], :] += wfv1 * N[i, :]
		VertexNormals[FV.faces[i][1], :] += wfv2 * N[i, :]
		VertexNormals[FV.faces[i][2], :] += wfv3 * N[i, :]
		
		"""
		Calculate areas for weights according to Mayar et al. [2002]
		Check if the triangle is obtuse, right or acute
		"""
		"Changed shape for ew"
		
		if ew[0][i] <= 0:
			Acorner[i][1] = -0.25 * l2[i][2] * Af[i] / (np.dot(e0[i, :], np.transpose(e2[i, :])))
			Acorner[i][2] = -0.25 * l2[i][1] * Af[i] / (np.dot(e0[i, :], np.transpose(e1[i, :])))
			Acorner[i][0] = Af[i] - Acorner[i][2] - Acorner[i][1]
		elif ew[1][i] <= 0:
			Acorner[i][2] = -0.25 * l2[i][0] * Af[i] / (np.dot(e1[i, :], np.transpose(e0[i, :])))
			Acorner[i][0] = -0.25 * l2[i][2] * Af[i] / (np.dot(e1[i, :], np.transpose(e2[i, :])))
			Acorner[i][1] = Af[i] - Acorner[i][2] - Acorner[i][0]
		elif ew[2][i] <= 0:
			Acorner[i][0] = -0.25 * l2[i][1] * Af[i] / (np.dot(e2[i, :], np.transpose(e1[i, :])))
			Acorner[i][1] = -0.25 * l2[i][0] * Af[i] / (np.dot(e2[i, :], np.transpose(e0[i, :])))
			Acorner[i][2] = Af[i] - Acorner[i][1] - Acorner[i][0]
		else:
			ewscale = 0.5 * Af[i] / (ew[0][i] + ew[1][i] + ew[2][i])
			Acorner[i][0] = ewscale * (ew[1][i] + ew[2][i])
			Acorner[i][1] = ewscale * (ew[0][i] + ew[2][i])
			Acorner[i][2] = ewscale * (ew[1][i] + ew[0][i])
		
		Avertex[FV.faces[i][0]] += Acorner[i][0]
		Avertex[FV.faces[i][1]] += Acorner[i][1]
		Avertex[FV.faces[i][2]] += Acorner[i][2]
		
		" Calcul initial coordinate system "
		up[FV.faces[i][0], :] = e2_norm[i, :]
		up[FV.faces[i][1], :] = e0_norm[i, :]
		up[FV.faces[i][2], :] = e1_norm[i, :]
	
	VertexNormals = normr(VertexNormals)
	
	" Calcul initial vertex coordinate system"
	
	for i in range(np.shape(FV.vertices)[0]):
		up[i, :] = np.cross(up[i, :], VertexNormals[i, :])
		up[i, :] = up[i, :] / np.linalg.norm(up[i, :])
		vp[i, :] = np.cross(VertexNormals[i, :], up[i, :])
	
	print("Finished calculating vertex normals")
	return VertexNormals, Avertex, Acorner, up, vp


def getPrincipalCurvatures(FV, VertexSFM, up, vp):
	print("Calculating Principal Components ... Please wait")
	
	"Calculate principal curvatures"
	PrincipalCurvature = np.zeros((2, np.shape(FV.vertices)[0]))
	PrincipalDi1, PrincipalDi2 = [np.zeros((np.shape(FV.vertices)[0], 3)), np.zeros((np.shape(FV.vertices)[0], 3))]
	for i in range(np.shape(FV.vertices)[0]):
		npp = np.cross(up[i, :], vp[i, :])
		r_old_u, r_old_v = RotateCoordinateSystem(up[i, :], vp[i, :], npp)
		ku = VertexSFM[i][0][0]
		kuv = VertexSFM[i][0][1]
		kv = VertexSFM[i][1][1]
		c, s, tt = 1, 0, 0
		if kuv != 0:
			"Jacobi rotation to diagonalize"
			h = 0.5 * (kv - ku) / kuv
			if h < 0:
				tt = 1 / (h - np.sqrt(1 + h ** 2))
			else:
				tt = 1 / (h + np.sqrt(1 + h ** 2))
			c = 1 / np.sqrt(1 + tt ** 2)
			s = tt * c
		k1 = ku - tt * kuv
		k2 = kv + tt * kuv
		if abs(k1) >= abs(k2):
			PrincipalDi1[i, :] = c * r_old_u - s * r_old_v
		else:
			[k1, k2] = [k2, k1]
			PrincipalDi1[i, :] = c * r_old_u + s * r_old_v
		PrincipalDi2[i, :] = np.cross(npp, PrincipalDi1[i, :])
		PrincipalCurvature[0][i] = k1
		PrincipalCurvature[1][i] = k2
		
		if np.isnan(k1) or np.isnan(k2):
			print("Nan")
	
	print("Finished Calculating principal components")
	return PrincipalCurvature, PrincipalDi1, PrincipalDi2


def RotateCoordinateSystem(up, vp, nf):
	r_new_u = up
	r_new_v = vp
	npp = np.cross(up, vp) / np.linalg.norm(np.cross(up, vp))
	ndot = np.dot(nf, np.transpose(npp))
	if ndot <= -1:
		r_new_u = -r_new_u
		r_new_v = -r_new_v
	perp = nf - ndot * npp
	dperp = (npp + nf) / (1 + ndot)
	r_new_u = r_new_u - dperp * np.dot(perp, np.transpose(r_new_u))
	r_new_v = r_new_v - dperp * np.dot(perp, np.transpose(r_new_v))
	
	return r_new_u, r_new_v

def somme_colonnes(X):
	"""
	
	:param X: A matrix
	:return:  A row that contains the sum of each column
	"""
	xx = list()
	for i in range(np.shape(X)[1]):
		xx.append(sum(X[:, i]))
	return np.array(xx)





if __name__ == '__main__':
	"Ellipse example"
	"""
	mesh = sio.load_mesh('ellipsoide.gii')
	"print(mesh.vertices)"
	"print(mesh.faces)"
	
	PrincipalCurvatures, PrincipalDir1, PrincipalDir2 = GetCurvaturesAndDerivatives(mesh)
	print("PrincipalCurvatures", np.shape(PrincipalCurvatures))
	print("Principale dir 1", np.shape(PrincipalDir1))
	print("Principale dir 2", np.shape(PrincipalDir2))
	gaussian_curv = PrincipalCurvatures[0, :] * PrincipalCurvatures[1, :]
	mean_curv = 0.5 * (PrincipalCurvatures[0, :] + PrincipalCurvatures[1, :])
	# {print("gaussian_curv", gaussian_curv)
	# print("mean_curv", mean_curv)
	# print("mesh.vertex_normals", mesh.vertex_normals)
	
	M = np.array([[1., 2., 3., 4.], [1., 2., 3., 4.]])
	x = np.transpose(np.sqrt(1 / somme_colonnes(np.transpose(M * M))))
	print(x)
	print(np.shape(x))
	print(np.reshape(x, (2, 1)))
	print(normr(M))
	mesh.show()
	pyglet_plot(mesh, np.transpose(gaussian_curv))
	pyglet_plot(mesh, mean_curv)
	
	"Sphere example"
	mesh = sio.load_mesh('sphere.gii')
	PrincipalCurvatures, PrincipalDir1, PrincipalDir2 = GetCurvaturesAndDerivatives(mesh)
	print("PrincipalCurvatures", np.shape(PrincipalCurvatures))
	print("Principale dir 1", np.shape(PrincipalDir1))
	print("Principale dir 2", np.shape(PrincipalDir2))
	gaussian_curv = PrincipalCurvatures[0, :] * PrincipalCurvatures[1, :]
	mean_curv = 0.5 * (PrincipalCurvatures[0, :] + PrincipalCurvatures[1, :])
	mesh.show()
	pyglet_plot(mesh, np.transpose(gaussian_curv))
	pyglet_plot(mesh, mean_curv)
	
	"Quad  100 example"
	mesh = sio.load_mesh('quadric_K1_1_K2_-1_10.gii')
	PrincipalCurvatures, PrincipalDir1, PrincipalDir2 = GetCurvaturesAndDerivatives(mesh)
	print("PrincipalCurvatures", np.shape(PrincipalCurvatures))
	print("Principale dir 1", np.shape(PrincipalDir1))
	print("Principale dir 2", np.shape(PrincipalDir2))
	gaussian_curv = PrincipalCurvatures[0, :] * PrincipalCurvatures[1, :]
	mean_curv = 0.5 * (PrincipalCurvatures[0, :] + PrincipalCurvatures[1, :])
	mesh.show()
	pyglet_plot(mesh, np.transpose(gaussian_curv))
	pyglet_plot(mesh, mean_curv)
	
	"Quad  2500 example"
	mesh = sio.load_mesh('quadric_K1_0_K2_1.gii')
	PrincipalCurvatures, PrincipalDir1, PrincipalDir2 = GetCurvaturesAndDerivatives(mesh)
	print("PrincipalCurvatures", np.shape(PrincipalCurvatures))
	print("Principale dir 1", np.shape(PrincipalDir1))
	print("Principale dir 2", np.shape(PrincipalDir2))
	gaussian_curv = PrincipalCurvatures[0, :] * PrincipalCurvatures[1, :]
	mean_curv = 0.5 * (PrincipalCurvatures[0, :] + PrincipalCurvatures[1, :])
	mesh.show()
	pyglet_plot(mesh, np.transpose(gaussian_curv))
	pyglet_plot(mesh, mean_curv)
	
	"Quad  225000 example"
	mesh = sio.load_mesh('quadric_K1_-1_K2_-1_150.gii')
	PrincipalCurvatures, PrincipalDir1, PrincipalDir2 = GetCurvaturesAndDerivatives(mesh)
	print("PrincipalCurvatures", np.shape(PrincipalCurvatures))
	print("Principale dir 1", np.shape(PrincipalDir1))
	print("Principale dir 2", np.shape(PrincipalDir2))
	gaussian_curv = PrincipalCurvatures[0, :] * PrincipalCurvatures[1, :]
	mean_curv = 0.5 * (PrincipalCurvatures[0, :] + PrincipalCurvatures[1, :])
	mesh.show()
	pyglet_plot(mesh, np.transpose(gaussian_curv))
	pyglet_plot(mesh, mean_curv)
	
	white_left_327680 sample
	mesh = sio.load_mesh('OAS1_0006_Lhemi.gii')
	PrincipalCurvatures, PrincipalDir1, PrincipalDir2 = GetCurvaturesAndDerivatives(mesh)
	gaussian_curv = PrincipalCurvatures[0, :] * PrincipalCurvatures[1, :]
	mean_curv = 0.5 * (PrincipalCurvatures[0, :] + PrincipalCurvatures[1, :])
	print(gaussian_curv)
	print(mean_curv)
	mesh.show()
	print("min cmean", min(mean_curv))
	print("max cmean", max(mean_curv))
	print("min cgauss", min(gaussian_curv))
	print("max cgauss", max(gaussian_curv))
	pyglet_plot(mesh, gaussian_curv)
	pyglet_plot(mesh, mean_curv)
"""
	t1 = time.time()
	fd = '/hpc/meca/users/bohi.a/Data/Models/week23_ref_surf/B0.txt'
	coords = np.loadtxt(fd, skiprows=1, max_rows=50943)
	faces = np.loadtxt(fd, skiprows=50945, dtype=np.int)
	mesh = trimesh.Trimesh(faces=faces - 1, vertices=coords, process=False)
	mesh.show()
	PrincipalCurvatures, PrincipalDir1, PrincipalDir2 = GetCurvaturesAndDerivatives(mesh)
	gaussian_curv = PrincipalCurvatures[0, :] * PrincipalCurvatures[1, :]
	mean_curv = 0.5 * (PrincipalCurvatures[0, :] + PrincipalCurvatures[1, :])
	print("le calcul a duré : ", time.time() - t1)
	print(gaussian_curv)
	print(mean_curv)
	mesh.show()
	print("min cmean", min(mean_curv))
	print("max cmean", max(mean_curv))
	print("min cgauss", min(gaussian_curv))
	print("max cgauss", max(gaussian_curv))
	vect_col_map = \
		trimesh.visual.color.interpolate(mean_curv, color_map='jet')
	
	if mean_curv.shape[0] == mesh.vertices.shape[0]:
		mesh.visual.vertex_colors = vect_col_map
	elif mean_curv.shape[0] == mesh.faces.shape[0]:
		mesh.visual.face_colors = vect_col_map
	mesh.show(background=[0, 0, 0, 255])
	
	mpl.rcParams.update({'font.size': 20})
	# fig = plt.figure(figsize=(8, 2))
	# ax1 = fig.add_axes([0.05, 0.80, 0.9, 0.15])
	fig, ax = plt.subplots(1, 1)
	ax.set_title("Colormap")
	
	vmin = np.min(mean_curv)
	vmax = np.max(mean_curv)
	norm = mpl.colors.Normalize(vmin=vmin, vmax=vmax)
	
	mpl.colorbar.ColorbarBase(ax, norm=norm, orientation='horizontal')
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
