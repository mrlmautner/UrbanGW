# -*- coding: utf-8 -*-
"""
Created on Mon Feb 18 20:49:58 2019

@author: MM
"""

import numpy as np

def nondom_sort(solutions):
    objs = solutions.copy()
    num_solutions = len(objs)

    # use a boolean index to keep track of nondominated solns
    keep = np.zeros(num_solutions, dtype = bool)

    for i in range(num_solutions):
        for j in range(num_solutions):
            a = objs[i]
            b = objs[j]
            if np.all(a <= b) & np.any(a < b):
                keep[i] = True
                keep[j] = False

    return keep