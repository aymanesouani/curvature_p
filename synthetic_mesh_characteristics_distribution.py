import matplotlib
matplotlib.use('TkAgg')
import matplotlib.pyplot as plt
import os
import trimesh
import numpy as np
import trimesh_caracteristics as tc
from trimesh_extension import trimesh_create_parametric_surfaces as gp



def trimesh_quad(K, nstep=10):
    X, Y, faces, Zs = gp.generate_quadric([K], nstep)
    Z =Zs[0]
    coords = np.array([X,Y,Z]).transpose()

    return trimesh.Trimesh(faces=faces, vertices=coords, process=False)

if __name__ == '__main__':
    nstep=20
    output_folder = '/hpc/meca/users/auzias/test_curvature'
    file_fig = os.path.join(output_folder,'synthethic_'+str(nstep)+'mesh_distrib_all.png')
    file_fig2 = os.path.join(output_folder,'synthetic_'+str(nstep)+'mesh_distrib.png')

    nb_quadrics = 100
    K = [1, 1]
    pop_mesh_carac = list()
    for quad_num in range(nb_quadrics):
        mesh = trimesh_quad(K,nstep)
        #mesh.show()
        pop_mesh_carac.append(tc.get_mesh_chracteristics(mesh))

    tc.plot_ditributions(np.array(pop_mesh_carac), file_fig, file_fig2)

    plt.show()