import sys
import os
import slam.io as sio
import slam.remeshing as srem
from slam import texture as stex
#import slam.plot as splt

if __name__ == '__main__':
    output_folder = '/hpc/meca/users/auzias/test_curvature/real_mesh_all_curvatures'
    target_spherical_mesh = sio.load_mesh('/hpc/meca/softs/dev/auzias/FSsphere.ico7.gii')
    if len(sys.argv) > 1:
        db_path = sys.argv[1]
        curv_types = [sys.argv[2]]
    else:
        #db_path = '/hpc/scalp/data/REPRO_database/FS_database_KKI_test_retest_FS6.0'
        db_path = '/hpc/scalp/data/REPRO_database/FS_database_OASIS_test_retest_FS5.3.0'
        curv_types = ['mean_Rusinkiewicz', 'mean_Dong', 'mean_Maillot', 'mean_PetitJean', 'mean_Peyre', 'mean_Taubin']#'gauss_Dong', 'gauss_Maillot', 'gauss_PetitJean', 'gauss_Peyre', 'gauss_Rusinkiewicz', 'gauss_Taubin']



    subj_files_list=os.listdir(db_path)
    subjects_list = list()
    subjects_mesh_files = list()
    for fil in subj_files_list:
        if fil.find('.') == -1 and fil.find('sanlm') == -1:
            lh_mesh_file = os.path.join(db_path, fil, 'surf', 'lh.sphere.reg.gii')
            rh_mesh_file = os.path.join(db_path, fil, 'surf', 'rh.sphere.reg.gii')
            if os.path.exists(lh_mesh_file) and os.path.exists(rh_mesh_file):
                subjects_list.append(fil)
                subjects_mesh_files.append(rh_mesh_file)
                subjects_mesh_files.append(lh_mesh_file)
    print('nb of meshes found: '+str(len(subjects_mesh_files)))

    for curv_type in curv_types:
        mesh_files = list()
        tex_files = list()
        for subject_mesh_file in subjects_mesh_files:
            filename_split = subject_mesh_file.split('/')
            subj = filename_split[6]
            filename = filename_split[-1]
            side = filename[0]
            print(subj+' '+side)
            curv = os.path.join(output_folder, subj + '.'+side+'h.white_curv_'+curv_type+'_NN_FSsphere.ico7.gii')
            if not os.path.exists(curv):
                mesh_files.append(subject_mesh_file)
                tex_files.append(curv)

        print('nb of meshes to be processed: '+str(len(mesh_files)))
        #mesh_files_i = mesh_files[:2]
        #tex_files_i = tex_files[:2]
        for mesh_file, tex_file in zip(mesh_files, tex_files):
            print('remeshing for '+tex_file)
            file_tex_in = tex_file[:-21]+'.gii'
            print(file_tex_in)

            source_tex = sio.load_texture(file_tex_in)
            source_spherical_mesh = sio.load_mesh(mesh_file)
            interpolated_tex = \
                srem.spherical_interpolation_nearest_neigbhor(source_spherical_mesh, target_spherical_mesh, source_tex.darray)
            #splt.pyglet_plot(target_spherical_mesh, interpolated_tex, -1, 1)

            output_texture = stex.TextureND(darray=interpolated_tex)
            sio.write_texture(output_texture, tex_file)
