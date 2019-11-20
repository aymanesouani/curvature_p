import matplotlib
matplotlib.use('TkAgg')
import matplotlib.pyplot as plt
import os
import trimesh_caracteristics as tc
# import trimesh_gifti as tg
import numpy as np
import  slam as sio

if __name__ == '__main__':
    output_folder = '/hpc/meca/users/auzias/test_curvature'
    db_path = '/hpc/meca/data/OASIS/BV_OASIS/OASIS'

    file_fig = os.path.join(output_folder,'OASIS_BV_mesh_distrib_all.png')
    file_fig2 = os.path.join(output_folder,'OASIS_BV_mesh_distrib.png')
    FS_mesh_path = 't1mri/freesurfer/FS_import_BV_mesh/segmentation/mesh'

    #file_fig = os.path.join(output_folder,'OASIS_FS_mesh_distrib_all.png')
    #file_fig2 = os.path.join(output_folder,'OASIS_FS_mesh_distrib.png')
    #FS_mesh_path = 't1mri/freesurfer/FS_import_FS_mesh/segmentation/mesh'


    BV_sides=['L','R']
    subjects_list1 = list()
    subj_files_list=os.listdir(db_path)
    for fil in subj_files_list:
        if fil.find('.') == -1:
            subjects_list1.append(fil)
    print('nb of subjects to be processed : '+str(len(subjects_list1)))
    print(subjects_list1)
    subjects_list = subjects_list1

    side = BV_sides[0]
    mesh_path = FS_mesh_path
    pop_mesh_carac = list()
    for subject in subjects_list:
        print(subject)
        mesh_file = os.path.join(db_path, subject, mesh_path, subject+'_'+side+'white.gii')
        mesh = sio.load(mesh_file)
        pop_mesh_carac.append(tc.get_mesh_chracteristics(mesh))

    tc.plot_ditributions(np.array(pop_mesh_carac), file_fig, file_fig2)

    plt.show()