# -*- coding: utf-8 -*-
"""
Created on Mon Jan  6 13:05:55 2020

@author: MM
"""
import plot_results as pltvm
#df_Bdget, mthly_Bdget, df_CumSum = pltvm.get_budgets(scenario_names, mapTitles, s_heads)
pltvm.plt_cum_sum('Combined_20200106', scenario_names, mapTitles, df_CumSum, start='01-31-1985', end='12-31-2013')