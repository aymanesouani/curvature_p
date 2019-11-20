# -*- coding: utf-8 -*-
"""
Created on Tue Mar 19 13:53:29 2013

@author: auzias
"""

import matplotlib.pyplot as plt
import os
import numpy as np
import Gradient as gd
from soma import aims

if __name__ == "__main__":

    fs_database_path ='/hpc/scalp/data/REPRO_database/FS_database_KKI_test_retest_FS6.0'
    db_pits_path = '/hpc/scalp/data/REPRO_database/SulcalPits'
    BV_sides=['L','R']
    sides=['lh','rh']
    subjects_list = list()
    subj_files_list=os.listdir(db_pits_path)
    for fil in subj_files_list:
        if fil.find('.') == -1:
            subjects_list.append(fil)
    print('nb of subjects to be processed : '+str(len(subjects_list)))
    print(subjects_list)

    side = sides[0]
    BV_side = BV_sides[0]

    mesh_path = 'surf'
    mesh_name = '.white.gii'


    subj = subjects_list[0]

    DPF_file = os.path.join(db_pits_path,subj,subj+'_'+BV_side+'white_DPF.gii')
    mesh_file = os.path.join(fs_database_path, subj, mesh_path, side+mesh_name)
    tex_DPF = aims.read(DPF_file)
    a_tex_DPF = np.array(tex_DPF[0])

    Grad = gd.Gradient(aims.read(mesh_file), aims.read(DPF_file))
    vectGrad=np.array(Grad.values())

    normVectGrad = np.sqrt(np.sum(np.power(vectGrad,2), axis=1))
    normVectGrad.shape

    tex_out = aims.TimeTexture_FLOAT(1, normVectGrad.shape[0])
    tex_out[0].assign(normVectGrad)
    aims.write(tex_out, os.path.join(db_pits_path,subj,subj+'_'+BV_side+'white_DPF_gradient.gii'))

    f , ax = plt.subplots(1,1)
    ax.hist(a_tex_DPF,100)
    plt.show()

    f , ax = plt.subplots(1,1)
    ax.hist(normVectGrad,100)
    plt.show()

    f, ax = plt.subplots(1,1)
    ax.plot(a_tex_DPF, normVectGrad,'o')
    plt.show()

    y=np.zeros(normVectGrad.shape)
    y[a_tex_DPF>0]=1
    f, ax = plt.subplots(1,1)
    ax.plot(y, normVectGrad,'o')
    plt.show()

    DPF_thresh = np.median(a_tex_DPF)
    dat = [normVectGrad[a_tex_DPF>DPF_thresh],normVectGrad[a_tex_DPF<DPF_thresh]]
    f, ax = plt.subplots(1,1)
    ax.boxplot(dat)
    plt.show()


    OASIS_diff = []
    OASIS_subjects_list = []
    KKI_diff = []
    KKI_subjects_list = []
    for ind_s,s in enumerate(subjects_list):
        if 'MR2' in s:
            if 'OAS1' in s:
                if s[:-4] in subjects_list:
                    corresp_ind = subjects_list.index(s[:-4])
                    corresp_s = subjects_list[corresp_ind]
                    MR2_tex_DPF = gi.read(os.path.join(db_pits_path,s,s+'_'+side+'white_DPF_FSsphere.ico7.gii'))
                    MR1_tex_DPF = gi.read(os.path.join(db_pits_path,corresp_s,corresp_s+'_'+side+'white_DPF_FSsphere.ico7.gii'))
                    OASIS_diff.append(MR2_tex_DPF.darrays[0].data - MR1_tex_DPF.darrays[0].data)
                    OASIS_subjects_list.append(s[:-4])
            else:
                if s[:-4]+'_MR1' in subjects_list:
                    corresp_ind = subjects_list.index(s[:-4]+'_MR1')
                    corresp_s = subjects_list[corresp_ind]
                    MR2_tex_DPF = gi.read(os.path.join(db_pits_path,s,s+'_'+side+'white_DPF_FSsphere.ico7.gii'))
                    MR1_tex_DPF = gi.read(os.path.join(db_pits_path,corresp_s,corresp_s+'_'+side+'white_DPF_FSsphere.ico7.gii'))
                    KKI_diff.append(MR2_tex_DPF.darrays[0].data - MR1_tex_DPF.darrays[0].data)
                    KKI_subjects_list.append(s[:-4])



    KKI_diff = np.array(KKI_diff).squeeze()
    KKI_sum_diff = np.sum(np.abs(KKI_diff),1)/KKI_diff.shape[1]

    OASIS_diff = np.array(OASIS_diff).squeeze()
    OASIS_sum_diff = np.sum(np.abs(OASIS_diff),1)/OASIS_diff.shape[1]
    print(np.abs(OASIS_diff).shape)
    print(OASIS_sum_diff.shape)
    f, ax = plt.subplots(1,1)
    f.suptitle('mean(abs(DPF_MR1-DPF_MR2))')
    x = np.ones(len(OASIS_sum_diff))
    ax.plot(x,OASIS_sum_diff,'ob')
    x = 2*np.ones(len(KKI_sum_diff))
    ax.plot(x,KKI_sum_diff,'or')
    ax.grid(True)
    ax.set_xlim((0, 3))
    ax.legend(['OASIS','KKI'])
    f.set_size_inches(10.5, 10.5)
    plt.savefig(file_fig, bbox_inches='tight')#, dpi=300)



    f, ax = plt.subplots(1,1)
    f.suptitle('KKI')
    ax.hist(KKI_diff.flatten(),100)
    f, ax = plt.subplots(1,1)
    f.suptitle('OASIS')
    ax.hist(OASIS_diff.flatten(),100)

    plt.show()
