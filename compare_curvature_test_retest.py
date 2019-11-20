import os
import slam.io as sio
from scipy.stats import pearsonr
import numpy as np
import matplotlib
matplotlib.use('TkAgg')
import matplotlib.pyplot as plt

if __name__ == '__main__':
    output_folder = '/hpc/meca/users/auzias/test_curvature/real_mesh_all_curvatures'
    target_spherical_mesh = sio.load_mesh('/hpc/meca/softs/dev/auzias/FSsphere.ico7.gii')

    curv_types = ['mean_Rusinkiewicz', 'mean_Dong', 'mean_PetitJean', 'mean_Peyre', 'mean_Taubin', 'gauss_Dong', 'gauss_PetitJean', 'gauss_Peyre', 'gauss_Rusinkiewicz', 'gauss_Taubin']
    #'mean_Maillot', 'gauss_Maillot'
    sides = ['lh']

    files_list=os.listdir(output_folder)
    correl = list()
    diff = list()
    for curv_type in curv_types:
        print(curv_type+'------------------------------------------------')
        for side in sides:
            tex_files = list()
            subjects = list()
            for fil in files_list:
                if curv_type in fil and 'NN_FSsphere.ico7.gii' in fil and side in fil:
                    tex_files.append(fil)
                    filename_split = fil.split('.')
                    subj = filename_split[0]
                    filename = filename_split[-1]
                    side = filename_split[1]
                    subjects.append(subj)

            #print('nb of tex found: ' + str(len(tex_files)))
            #print(subjects)

            OASIS_diff = []
            OASIS_corr_coef = []
            OASIS_subjects_list = []
            KKI_diff = []
            KKI_corr_coef = []
            KKI_subjects_list = []
            for ind_s,s in enumerate(subjects):
                if 'MR2' in s:
                    if 'OAS1' in s:
                        if s[:-4] in subjects:
                            corresp_ind = subjects.index(s[:-4])
                            corresp_s = subjects[corresp_ind]
                            MR2_tex = sio.load_texture(os.path.join(output_folder, tex_files[ind_s]))
                            MR1_tex = sio.load_texture(os.path.join(output_folder, tex_files[corresp_ind]))
                            OASIS_diff.append(MR2_tex.darray - MR1_tex.darray)
                            OASIS_subjects_list.append(s[:-4])
                            per = pearsonr(MR2_tex.darray, MR1_tex.darray)
                            OASIS_corr_coef.append(per[0][0])
                    else:
                        if s[:-4]+'_MR1' in subjects:
                            corresp_ind = subjects.index(s[:-4]+'_MR1')
                            corresp_s = subjects[corresp_ind]
                            MR2_tex = sio.load_texture(os.path.join(output_folder, tex_files[ind_s]))
                            MR1_tex = sio.load_texture(os.path.join(output_folder, tex_files[corresp_ind]))
                            KKI_diff.append(MR2_tex.darray - MR1_tex.darray)
                            KKI_subjects_list.append(s[:-4])
                            per = pearsonr(MR2_tex.darray, MR1_tex.darray)
                            KKI_corr_coef.append(per[0][0])


            correl.append([np.mean(OASIS_corr_coef), np.mean(KKI_corr_coef)])
            diff.append([np.mean(np.abs(OASIS_diff)), np.mean(np.abs(KKI_diff))])

    for curv_type, curv_corr, curv_diff in zip(curv_types, correl, diff):
        print(curv_type )
        print('average correl OASIS : '+str(curv_corr[0]))
        print('average correl KKI : '+str(curv_corr[1]))
        print('average diff OASIS : '+str(curv_diff[0]))
        print('average diff KKI : '+str(curv_diff[1]))


    plt.figure()
    plt.plot(MR2_tex.darray, MR1_tex.darray,'.')
    plt.figure()
    plt.hist(KKI_diff[-1], 100)
    plt.show()