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
    mat_file.write("file_out='%s';\n" %file_out)
    mat_file.write("folder_out='%s';\n" % folder_out)
    mat_file.write("compute_all_curvatures(file_mesh_in, file_out, folder_out)\n")
    mat_file.write("exit\n")
    mat_file.close()

    cmd = matlab_exe + " -nodisplay -r "+matfilePath#+"'"
    print(cmd)
    os.system(cmd)
    os.chdir(script_path)
    # print(os.getcwd())
    # os.popen(cmd)
    # p = subprocess.Popen([matlab_exe, '-r', matfilePath])
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

"""def caret_curvatures(mesh_file, file_out_name):
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
    output_folder = '/hpc/meca/users/auzias/test_curvature'

    mesh_files = []
    Ks = [[1, 1], [1, 0], [0, 0], [-1, 1]]
    for K in Ks:
        mesh_files.append(os.path.join(output_folder,'quadric_K1_'+str(K[0])+'_K2_'+str(K[1])+'.gii'))
    # mesh_files.append(os.path.join(output_folder,'FSsphere.ico7.gii'))
    # mesh_files.append(os.path.join(output_folder,'KKI2009_113_MR1.lh.white.gii'))
    for mesh_file in mesh_files:
        print('computing curvatures for '+mesh_file)
        file_out_name = mesh_file[:-4]+'_curv'
        # print('--caret')
        # caret_curvatures(mesh_file, file_out_name)
        # print('--brainvisa')
        # brainvisa_curvatures(mesh_file, file_out_name)
        print('--matlab')
        matlab_curvatures(mesh_file, output_folder)
     """