import os
import sys
import subprocess
import time

def matlab_curvatures(mesh_file, file_out, folder_out):

    matlab_exe = "/hpc/soft/matlab/matlab2014a/bin/matlab" #"/usr/bin/matlab-R2013a"#


    script_path = '/hpc/meca/softs/dev/auzias/matlab/surface_process'
    c=time.time()
    matfilePath = 'tmp_compute_curvature'+str(int(1000000*c))
    print(matfilePath)
    mat_file = open(matfilePath+'.m', 'w')
    #mat_file.write("cd('%s');\n" % script_path)
    #mat_file.write("pwd\n")
    mat_file.write("addpath(genpath('%s'));\n" % script_path)
    #mat_file.write("addpath('io');\n")
    mat_file.write("file_mesh_in='%s';\n" % mesh_file)
    mat_file.write("folder_out='%s';\n" % folder_out)
    mat_file.write("file_out='%s';\n" % file_out)
    mat_file.write("compute_all_curvatures(file_mesh_in,file_out,folder_out)\n")
    mat_file.write("exit\n")
    mat_file.close()

    cmd = matlab_exe + " -nodisplay -r "+matfilePath#+"'"
    print(cmd)
    os.system(cmd)
    os.chdir(script_path)
    #print(os.getcwd())
    #os.popen(cmd)
    #p = sub.Popen([matlab_exe, '-r', matfilePath])
    # a = subprocess.getoutput(cmd)
    # print(a)
    # os.remove(matfilePath+'.m')

def brainvisa_curvatures(mesh_file, file_out_name):
    bv_cmd = '/hpc/meca/users/auzias/brainvisa-4.5_with_patch/bin/AimsMeshCurvature'
    # FEM method
    file_out = file_out_name+'_mean_fem.gii'
    cmd = bv_cmd+' -i '+mesh_file+' -o '+file_out+' -m fem -r 1'
    a = subprocess.getoutput(cmd)
    print(a)
    # Boix method
    file_out = file_out_name+'_mean_boix.gii'
    cmd = bv_cmd+' -i '+mesh_file+' -o '+file_out+' -m boix -r 1'
    a = subprocess.getoutput(cmd)
    print(a)
    # Barycenter method
    file_out = file_out_name+'_mean_bary.gii'
    cmd = bv_cmd+' -i '+mesh_file+' -o '+file_out+' -m barycenter -r 1'
    a = subprocess.getoutput(cmd)
    print(a)
    # Boix Gaussian method
    file_out = file_out_name+'_gauss_boixGauss.gii'
    cmd = bv_cmd+' -i '+mesh_file+' -o '+file_out+' -m boixgaussian -r 1'
    a = subprocess.getoutput(cmd)
    print(a)

def caret_curvatures(mesh_file, file_out_name):
    workbench_cmd = '/hpc/soft/workbench/bin_linux64/'
    cmd = workbench_cmd+'wb_command -surface-curvature '+mesh_file+' -mean '+file_out_name+'_mean_Maillot.gii -gauss '+file_out_name+'_gauss_Maillot.gii'
    a = subprocess.getoutput(cmd)
    print(a)

def freesurfer_curvatures(mesh_file, file_out_name):
    fs_cmd = 'mris_curvature -w -max -min'
    cmd = fs_cmd+' '+mesh_file
    a = subprocess.getoutput(cmd)
    print(a)
    #convert to gifti
    #mris_convert -c /hpc/meca/users/auzias/test_curvature/unknown.quadric_K1_-1_K2_-1.gii.min /hpc/meca/users/auzias/test_curvature/quadric_K1_-1_K2_-1.gii  /hpc/meca/users/auzias/test_curvature/quadric_K1_-1_K2_-1.min.gii



if __name__ == '__main__':
    output_folder = '/hpc/meca/users/auzias/test_curvature/real_mesh_all_curvatures'
    #db_path = '/hpc/scalp/data/REPRO_database/FS_database_KKI_test_retest_FS6.0'
    db_path = '/hpc/scalp/data/REPRO_database/FS_database_OASIS_test_retest_FS5.3.0'

    subj_files_list=os.listdir(db_path)
    subjects_list = list()
    mesh_files = list()
    for fil in subj_files_list:
        if fil.find('.') == -1 and fil.find('sanlm') == -1:
            lh_mesh_file = os.path.join(db_path,fil,'surf','lh.white.gii')
            rh_mesh_file = os.path.join(db_path,fil,'surf','rh.white.gii')
            if os.path.exists(lh_mesh_file) and os.path.exists(rh_mesh_file):
                subjects_list.append(fil)
                last_curv_r = os.path.join(output_folder, fil+'.rh.white_curv_gauss_PetitJean.gii')
                if not os.path.exists(last_curv_r):
                    mesh_files.append(rh_mesh_file)
                last_curv_l = os.path.join(output_folder, fil+'.lh.white_curv_gauss_PetitJean.gii')
                if not os.path.exists(last_curv_l):
                    mesh_files.append(lh_mesh_file)


    print('nb of meshes to be processed : '+str(len(mesh_files)))
    print(mesh_files)
    for mesh_file in mesh_files:
        print('computing curvatures for '+mesh_file)
        filename_split = mesh_file.split('/')
        subj = filename_split[6]
        filename = filename_split[-1]
        print(subj)
        file_out_name = os.path.join(output_folder, subj+'.'+filename[:-4]+'_curv')
        print(file_out_name)
        print('--caret')
        #caret_curvatures(mesh_file, file_out_name)
        # print('--brainvisa')
        # brainvisa_curvatures(mesh_file, file_out_name)
        print('--matlab')
        file_out = subj+'.'+filename[:-4]
        #matlab_curvatures(mesh_file, file_out, output_folder)

