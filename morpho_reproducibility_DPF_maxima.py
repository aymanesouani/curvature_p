# -*- coding: utf-8 -*-
"""
Created on Tue Mar 19 13:53:29 2013

@author: auzias
"""

import matplotlib
matplotlib.use('TkAgg')
import matplotlib.pyplot as plt
import os
import numpy as np
from nibabel import gifti as gi
from ICC import ICC_rep_anova

def compute_ICC(KKI_data_test, KKI_data_retest):
    (KKI_nb_subj, KKI_nb_meas) = KKI_data_test.shape
    print('nb_subj='+str(KKI_nb_subj))
    print('nb_measures='+str(KKI_nb_meas))
    KKI_icc = np.zeros(KKI_nb_meas)
    KKI_session_F = np.zeros(KKI_nb_meas)
    KKI_session_var = np.zeros(KKI_nb_meas)
    KKI_subject_var = np.zeros(KKI_nb_meas)
    KKI_R2 = np.zeros(KKI_nb_meas)
    for st in range(KKI_nb_meas):
        Y = np.array([KKI_data_test[:, st], KKI_data_retest[:, st]]).transpose()
        if len(Y)>0:
            if np.sum(Y)>0:
                KKI_icc[st], KKI_subject_var[st], KKI_session_var[st], KKI_session_F[st], _, _ = ICC_rep_anova(Y)
                t = np.corrcoef(Y.transpose())
                KKI_R2[st] = t[0, 1] * t[0, 1]
            else:
                KKI_icc[st] = 2
                KKI_subject_var[st] = 0
                KKI_session_var[st] = 0
                KKI_session_F[st] = 0
        else:
            KKI_icc[st] = 0
            KKI_R2[st] = 0
    return KKI_icc, KKI_R2, KKI_session_F, KKI_session_var, KKI_subject_var


def plot_axis(ax, OASIS_data_test,OASIS_data_retest, xmin, xmax):
    (OASIS_nb_subj, OASIS_nb_meas) = OASIS_data_test.shape
    y = np.arange(OASIS_nb_meas) + 1
    for s in range(OASIS_nb_subj):
        for i in range(len(y)):
            ax.plot(OASIS_data_test[s, i], y[i] + 0.2, '.b')
            ax.plot(OASIS_data_retest[s, i], y[i] - 0.2, '.r')
            x_t =[OASIS_data_test[s, i], OASIS_data_retest[s, i]]
            y_t = [y[i] + 0.2, y[i] - 0.2]
            ax.plot(x_t, y_t, 'k')
    ax.grid(True)
    ax.set_yticks(y)
    #ax.set_yticklabels(['left', 'right'])
    ax.legend(['test', 'retest'])
    ax.set_xlim((xmin, xmax))
    return ax


def plot_axis_ICC(ax, OASIS_icc, OASIS_R2):
    # axes[1].boxplot(data_icc,vert=0)
    y = np.arange(len(OASIS_icc)) + 1
    ax.plot(OASIS_icc, y, 'go-')
    ax.plot(OASIS_R2, y, 'mo-')

    ICC_m = np.mean(OASIS_icc)
    ax.plot([ICC_m,ICC_m], [y[0],y[-1]], 'g--')
    R2_m = np.mean(OASIS_R2)
    ax.plot([R2_m,R2_m], [y[0],y[-1]], 'm--')

    ax.grid(True)
    ax.set_xlim((0, 1))


if __name__ == "__main__":


    db_pits_path = '/hpc/scalp/data/REPRO_database/SulcalPits'

    BV_sides=['L','R']
    subjects_list = list()
    subj_files_list=os.listdir(db_pits_path)
    for fil in subj_files_list:
        if fil.find('.') == -1:
            subjects_list.append(fil)
    print('nb of subjects to be processed : '+str(len(subjects_list)))
    print(subjects_list)

    side = BV_sides[0]

    file_fig_diff = '/hpc/scalp/data/REPRO_database/morphometry_DPF/'+side+'DPF_maxima_diff.png'
    file_fig_ICC = '/hpc/scalp/data/REPRO_database/morphometry_DPF/'+side+'DPF_maxima_ICC.png'

    OASIS_count = []
    OASIS_subjects_list = []
    KKI_count = []
    KKI_subjects_list = []
    for ind_s,s in enumerate(subjects_list):
        if 'MR2' in s:
            if 'OAS1' in s:
                if s[:-4] in subjects_list:
                    corresp_ind = subjects_list.index(s[:-4])
                    corresp_s = subjects_list[corresp_ind]
                    MR2_tex_DPF = gi.read(os.path.join(db_pits_path,s,s+'_'+side+'white_DPF_maxima_FSsphere.ico7.gii.gii'))
                    MR1_tex_DPF = gi.read(os.path.join(db_pits_path,corresp_s,corresp_s+'_'+side+'white_DPF_maxima_FSsphere.ico7.gii.gii'))
                    OASIS_count.append([np.sum(MR1_tex_DPF.darrays[0].data>0), np.sum(MR2_tex_DPF.darrays[0].data>0)])
                    OASIS_subjects_list.append(s[:-4])
            else:
                if s[:-4]+'_MR1' in subjects_list:
                    corresp_ind = subjects_list.index(s[:-4]+'_MR1')
                    corresp_s = subjects_list[corresp_ind]
                    MR2_tex_DPF = gi.read(os.path.join(db_pits_path,s,s+'_'+side+'white_DPF_maxima_FSsphere.ico7.gii.gii'))
                    MR1_tex_DPF = gi.read(os.path.join(db_pits_path,corresp_s,corresp_s+'_'+side+'white_DPF_maxima_FSsphere.ico7.gii.gii'))
                    KKI_count.append([np.sum(MR1_tex_DPF.darrays[0].data>0), np.sum(MR2_tex_DPF.darrays[0].data>0)])
                    KKI_subjects_list.append(s[:-4])

    KKI_count = np.array(KKI_count).squeeze()
    KKI_diff = KKI_count[:,0]-KKI_count[:,1]

    OASIS_count = np.array(OASIS_count).squeeze()
    OASIS_diff = OASIS_count[:,0]-OASIS_count[:,1]

    KKI_icc, KKI_R2, KKI_session_F, KKI_session_var, KKI_subject_var = compute_ICC(np.array([KKI_count[:,0]]).T, np.array([KKI_count[:,1]]).T)
    OASIS_icc, OASIS_R2, OASIS_session_F, OASIS_session_var, OASIS_subject_var = compute_ICC(np.array([OASIS_count[:,0]]).T, np.array([OASIS_count[:,1]]).T)


    xmax = np.max([np.max(KKI_count), np.max(OASIS_count)])+10
    xmin = np.min([np.min(KKI_count), np.min(OASIS_count)])-10

    f, axes = plt.subplots(2, 2, sharey=True)
    f.suptitle(side+' number of maxima in the DPF')
    plot_axis(axes[0, 0], np.array([OASIS_count[:,0]]).T, np.array([OASIS_count[:,1]]).T, xmin, xmax)
    axes[0, 0].set_title('OASIS '+side)
    plot_axis(axes[1, 0], np.array([KKI_count[:,0]]).T, np.array([KKI_count[:,1]]).T, xmin, xmax)
    axes[1, 0].set_title('KKI '+side)
    plot_axis_ICC(axes[0, 1], OASIS_icc, OASIS_R2)
    #axes[0, 1].set_title('OASIS')
    plot_axis_ICC(axes[1, 1], KKI_icc, KKI_R2)
    #axes[1, 1].set_title('KKI')
    f.set_size_inches(10.5, 10.5)
    plt.savefig(file_fig_ICC, bbox_inches='tight')#, dpi=300)


    print(OASIS_diff.shape)
    f, ax = plt.subplots(1,1)
    f.suptitle(side+' MR1-MR2 for the number of maxima in the DPF')
    x = np.ones(len(OASIS_count[:,0]))
    ax.plot(x,OASIS_diff,'ob')
    x = 2*np.ones(len(KKI_count[:,0]))
    ax.plot(x,KKI_diff,'or')
    ax.grid(True)
    ax.set_xlim((0, 3))
    ax.legend(['OASIS','KKI'])
    f.set_size_inches(10.5, 10.5)
    plt.savefig(file_fig_diff, bbox_inches='tight')#, dpi=300)
    plt.show()


