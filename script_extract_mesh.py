import sys
import os
import slam.io as sio
import nibabel.freesurfer.io as nbio
import numpy as np
from slam import topology as stop
import trimesh
import slam.plot as splt
import slam.remeshing as srem
from slam import texture as stex
import matplotlib.pyplot as plt
import scipy.sparse as sp
# from quadrics import *

# Data
data_path = '/hpc/meca/users/souani.a/data'
mesh_file = os.path.join(data_path,'lhwhite.gii')
mesh = sio.load_mesh(mesh_file)
texture_file = os.path.join(data_path,'OASIS_0004_Lwhite_bassins.gii')
texture = sio.load_texture(texture_file)
texture=texture.darray

# Mesh + texture
splt.pyglet_plot(mesh, texture, 'jet', True)

unique_values=np.unique(texture)
print(unique_values)
unique_values=np.setdiff1d(unique_values,0)
print(unique_values)
# Density
for i,val in enumerate(unique_values):
    inds = texture ==val
    sub_meshes, sub_tex, sub_corresp = stop.cut_mesh(mesh, inds)
    print(sub_meshes[0].vertices.shape[0])
    print(sub_meshes[0].area_faces.sum())
    print(sub_meshes[0].vertices.shape[0]/sub_meshes[0].area_faces.sum())

# Sub_meshes

#plt.figure(1)
#plt.figure(2)
for i,val in enumerate(unique_values):g
    inds = texture ==val
    sub_meshes, sub_tex, sub_corresp = stop.cut_mesh(mesh, inds)
    print("----------------------------------------------------")
    print(sub_corresp)
    print(sub_corresp[0][0])
    print(sub_corresp[0][-1])
    print("----------------------------------------------------")
    print(sub_meshes[0].vertices.shape)
    print(sub_meshes[0].area_faces.sum())
    # Distributions
    plt.figure(1)
    plt.subplot(2,3,i+1)
    nbins=np.sqrt(len(sub_meshes[0].edges_unique_length))
    plt.hist(sub_meshes[0].edges_unique_length,int(nbins))
    plt.figure(2)
    plt.subplot(2,3,i+1)
    nbins=np.sqrt(len(sub_meshes[0].area_faces))
    plt.hist(sub_meshes[0].area_faces,int(nbins))
    # Patches
    A=stop.edges_to_adjacency_matrix(sub_meshes[0])
    L=sp.diags([A.sum(axis=0)],[0],A.shape)-A
    L=L.toarray()
    curvature=np.sum(np.dot(L, sub_meshes[0].vertices)*sub_meshes[0].vertex_normals,axis=1)
    splt.pyglet_plot(sub_meshes[0], curvature, 'jet', True)
    sio.write_mesh(sub_meshes[0], "sub_mesh"+str(i))

plt.figure(1)
plt.title('Distribution of edges length',pad=170)
plt.figure(2)
plt.title('Distribution of faces area',pad=170)


# Density

all_areas=[]
all_pts=[]
for i,val in enumerate(unique_values):
    inds = texture == val
    sub_meshes, sub_tex, sub_corresp = stop.cut_mesh(mesh, inds)
    sub_meshes[0].area_faces.sum()
    all_areas.append(sub_meshes[0].area_faces.sum())
    all_pts.append(len(sub_meshes[0].vertices))

print(all_areas)
print(all_pts)

# Fitting meshes by quadrics

#for i,val in enumerate(unique_values):
#    inds = texture ==val
#    sub_meshes, sub_tex, sub_corresp = stop.cut_mesh(mesh, inds)
#    coeffs=estimate_coeffs(sub_meshes[0])
#    print(coeffs)