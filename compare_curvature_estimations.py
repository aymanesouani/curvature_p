import matplotlib
matplotlib.use('TkAgg')
import matplotlib.pyplot as plt
import nibabel as nb
import numpy as np
import os
import slam.generate_parametric_surfaces as gp
import slam.io as sio

if __name__ == '__main__':
    """main_folder = '/hpc/meca/users/auzias/test_curvature'
    file_sphere_mesh = os.path.join(main_folder,'FSsphere.ico7.gii')
    sphere_mesh = ng.read(file_sphere_mesh)
    mesh_coords = sphere_mesh.darrays[0].data
    bary = np.mean(mesh_coords)
    print(bary)
    print(np.min(mesh_coords))
    print(np.max(mesh_coords))"""

    main_folder = "/hpc/meca/users/souani.a/curvature/Compute_curvatures/Quadrics"
    curv_types = ['Peyre', 'Rusinkiewicz', 'PetitJean', 'Taubin', 'Dong']#, 'Peyre', 'fem', 'boix', 'bary']
    curv_types_leg = curv_types.copy()
    curv_types_leg.append('analytic')
    meshes = list()
    curv_mean = list()
    Ks = [[-1, 1], [1, 1], [0, 1]]
    for K in Ks:
        meshes.append('hex_quad_k1_'+str(K[0])+'_k2_'+str(K[1])+'_100')
        os.chdir(r"/hpc/meca/users/souani.a/curvature/Compute_curvatures/Quadrics")
        tmp_m = sio.load_mesh('hex_quad_k1_'+str(K[0])+'_k2_'+str(K[1])+'_100.gii')
        mesh_coords = tmp_m.vertices
        print(np.shape(mesh_coords))
        curv_mean.append(gp.quadric_curv_mean(K[0], K[1])(mesh_coords[:, 0], mesh_coords[:, 1]))
    #meshes.append('FSsphere.ico7')
    #meshes.append('KKI2009_113_MR1.lh.white')
    for mesh_type, curv_ref in zip(meshes, curv_mean):
        print(mesh_type)
        
        print("Shaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaape", np.shape(curv_ref))
        file_fig = os.path.join(main_folder, mesh_type+'_curv_mean_comparison_hist.png')
        curvatures = list()
        integral = list()
        for curv_type in curv_types:
            print('----'+curv_type)
            
            curv_tex_gii = nb.load(os.path.join("/hpc/meca/users/souani.a/curvature/Compute_curvatures", mesh_type+'_curv_mean_'+curv_type+'.gii'))

            curv_vals = curv_tex_gii.darrays[0].data.squeeze()
            if curv_type is 'PetitJean':
                curv_vals = curv_vals*2
            curvatures.append(curv_vals)
            print("Shaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaape :", np.shape(curv_vals))
            diff_curv = np.abs(curv_vals-curv_ref)
            print('  mean='+str(np.mean(curv_vals))+' std='+str(np.std(curv_vals)))
            print('  mean error='+str(np.mean(diff_curv))+' std error='+str(np.std(diff_curv)))
            integral.append(diff_curv)
            print('  integral of error = '+str(np.sum(diff_curv)) )
            print('  integral of normalized error = '+str(np.sum(diff_curv) / np.sum(np.abs(curv_ref))) )


        bins = np.linspace(-0.5, 0.5, 20)
        f, axes = plt.subplots(1, 1)
        f.suptitle('mean curvature '+mesh_type+' nbvert='+str(len(curv_vals)))
        for cc in curvatures:
            hist_t1, edges = np.histogram(cc, bins)
            left,right = edges[:-1],edges[1:]
            X = np.array([left,right]).T.flatten()
            Y = np.array([hist_t1,hist_t1]).T.flatten()
            axes.plot(X,Y)
        hist_ref, edges = np.histogram(curv_ref, bins)
        left,right = edges[:-1],edges[1:]
        X = np.array([left,right]).T.flatten()
        Y = np.array([hist_ref,hist_ref]).T.flatten()
        axes.plot(X,Y)

        axes.legend(curv_types_leg)
        axes.grid(True)

        f.set_size_inches(18.5, 10.5)
        plt.savefig(file_fig, bbox_inches='tight')#, dpi=300)

        bins = np.linspace(0, 0.1, 20)
        f, axes = plt.subplots(1, 1)
        f.suptitle('integral of error for curvature '+mesh_type+' nbvert='+str(len(curv_vals)))
        for cc in integral:
            hist_t1, edges = np.histogram(cc, bins)
            left,right = edges[:-1],edges[1:]
            X = np.array([left,right]).T.flatten()
            Y = np.array([hist_t1,hist_t1]).T.flatten()
            axes.plot(X,Y)

        axes.legend(curv_types)
        axes.grid(True)

    plt.show()
