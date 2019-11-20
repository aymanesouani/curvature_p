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

if __name__ == "__main__":

    BV_database_dir = '/hpc/scalp/data/REPRO_database/BV_database_KKI_test_retest/KKI2009/'
    acquisitions = ['default_acquisition', 'sanlm_denoised', 'FS_6_import_default_acquisition','FS_6_import_sanlm_denoised']
    sides = ['L','R']
    alphas = ['0.003','0.03','0.3']
    nb_subj=42
    bins = np.linspace(-10, 10, 100)

    subj_files_list = os.listdir(BV_database_dir)
    subj_list_processed = []
    for fil in subj_files_list:
        if not '.' in fil:
            subj_list_processed.append(fil[:-4])

    f1, axs = plt.subplots(len(alphas), 2, sharey=True, sharex=True)
    for indA, alpha in enumerate(alphas):

        for indS,side in enumerate(sides):
            list_DPF_diff_acq = []
            for acquisition in acquisitions:
                if acquisition is 'default_acquisition' or acquisition is 'sanlm_denoised':
                    analysis = 'morphologist_2015'
                else:
                    analysis = 'default_analysis'

                DPF_file_path = '/t1mri/'+acquisition+'/'+analysis+'/segmentation/mesh/surface_analysis/'

                # f1, axes = plt.subplots(1, nb_subj)
                # f1.suptitle('DPF for '+acquisition)
                list_subjects = []
                list_DPF_diff = []
                for ind,suj in enumerate(subj_list_processed[:nb_subj]):
                    file_path_MR1 = BV_database_dir + suj + '_MR1' + DPF_file_path + suj + '_MR1' + '_' + side + 'white_DPF_' + alpha + '.gii'
                    file_path_MR2 = BV_database_dir + suj + '_MR2' + DPF_file_path + suj + '_MR2' + '_' + side + 'white_DPF_' + alpha + '.gii'
                    #print(file_path_MR1)
                    if os.path.exists(file_path_MR1) and os.path.exists(file_path_MR2):
                        list_subjects.append(suj)
                        dpf_t1 = gi.read(file_path_MR1)
                        dpf_t2 = gi.read(file_path_MR2)
                        hist_t1, edges = np.histogram(dpf_t1.darrays[0].data, bins, normed=1)
                        # left,right = edges[:-1],edges[1:]
                        # X = np.array([left,right]).T.flatten()
                        # Y = np.array([hist_t1,hist_t1]).T.flatten()
                        # axes[ind].plot(X,Y,'b')
                        hist_t2, edges = np.histogram(dpf_t2.darrays[0].data, bins, normed=1)
                        # left,right = edges[:-1],edges[1:]
                        # X = np.array([left,right]).T.flatten()
                        # Y = np.array([hist_t2,hist_t2]).T.flatten()
                        # axes[ind].plot(X,Y,'r')
                        DPF_diff = np.sum(np.abs(hist_t1-hist_t2))
                        # axes[ind].legend(['t1','t2 sum_diff=%0.2f'+DPF_diff])

                        list_DPF_diff.append(DPF_diff)

                print('list_subjects : ', list_subjects)
                print('-------------------------------------')
                list_DPF_diff_acq.append(list_DPF_diff)
            list_DPF_diff_acq =np.array(list_DPF_diff_acq)

            axs[indA, indS].set_title('sum of bin-diff for side '+side+' for alpha='+alpha)
            axs[indA, indS].boxplot(list_DPF_diff_acq.T, vert=False)
            axs[indA, indS].set_yticklabels(acquisitions)
            axs[indA, indS].grid(True)

    plt.show()