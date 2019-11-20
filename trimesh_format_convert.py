import os
import slam.io as sio
import trimesh

if __name__ == '__main__':

    output_folder = '/hpc/meca/users/auzias/test_curvature'
    db_path = '/hpc/meca/data/OASIS/BV_OASIS/OASIS'
    FS_mesh_path = 't1mri/freesurfer/FS_import_FS_mesh/segmentation/mesh'


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
    subject = subjects_list[0]
    print(subject)
    mesh_file = os.path.join(db_path, subject, mesh_path, subject+'_'+side+'white.gii')
    mesh = sio.load(mesh_file)
    mesh.export(os.path.join(output_folder, subject+'_'+side+'white.ply'))
