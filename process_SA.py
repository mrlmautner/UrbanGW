# -*- coding: utf-8 -*-
"""
Created on Mon Apr 27 14:25:09 2020

@author: MM
"""

import plot_results as pltvm
import gwscripts.sensitivityanalysis.satools as sa 

# Set observation groups using k-means clustering based on observation characteristics (latitude, longitude, head elevation, and time [optional]). User chooses whether to use the absolute characteristic values (norm=False) or normalized based on the spread of the sample within each characteristic (norm=True)
df, obsinfo, obsstats, obsformation = pltvm.process_hobs('20200403_100000', '00000', obsinfo_loaded=True)
df.rename(columns={'absobserved': 'Z', 'time_series': 't'}, inplace=True)
df.set_index('obs_name', inplace=True)
group_label = sa.kmeans_obs(df, norm=True, time=False)
#group_label = sa.kmeans_obs(df)
#group_label = sa.kmeans_obs(df, time=False)
#group_label = sa.kmeans_obs(df, norm=True)

