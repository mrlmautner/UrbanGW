# -*- coding: utf-8 -*-
"""
Created on Mon Nov 16 17:07:37 2020

@author: MM
"""
import sa_analysis_functions as sa

sa_dataset = sa.dataset('20200921_100000-norm-5')
sa_dataset.plt_cluster()
sa_dataset.normalize_paramrange()
sa_dataset.cluster_kde(params=['HK_2','HK_3','HK_4','HK_5','SS_2','SY_2','SY_4','Q_3','LK_2','LK_3','TWU_1','TWU_2','TWU_3','RCH_2','IN_1'], threshold=0.05)
sa_dataset.objective_scatter(params=['HK_5'], threshold=0.05, show_all=True)
sa_dataset.historical_scatter(threshold=0.05, show_all=True)
sa_dataset.dominant_alternative_bar(threshold=0.05)
sa_dataset.param_sensitivity(threshold='0.05', constant='Alternative', def_alt=3)
sa_dataset.param_sensitivity(threshold='0.05', constant='Cluster', def_clust=5)
sa_dataset.param_sensitivity(threshold='0.05', constant='Objective', def_obj=0)
sa_dataset.obj_spread(threshold=0.05, max_val=80, step_val=5)
sa_dataset.heatmap_by_param(threshold=0.05, n_samples=100000, c='C-00005')