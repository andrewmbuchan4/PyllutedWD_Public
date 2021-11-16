#!/usr/bin/env python
# -*- coding: utf-8 -*-

from numba import jit

import numpy as np

@jit(nopython=True)
def frange(x, y, step):
    while x < y:
        yield x
        x += step

@jit(nopython=True)
def erfcc(x):
    z = abs(x)
    t = 1. / (1. + 0.5*z)
    r = t * np.exp(-z*z-1.26551223+t*(1.00002368+t*(.37409196+
        t*(.09678418+t*(-.18628806+t*(.27886807+
        t*(-1.13520398+t*(1.48851587+t*(-.82215223+
        t*.17087277)))))))))
    if (x >= 0.):
         return r
    else:
         return 2. - r

@jit(nopython=True)
def normcdf(x, mu, sigma):
    root_2 = 2.0**0.5
    t = x-mu;
    y = 0.5*erfcc(-t/(sigma*root_2));
    if y>1.0:
        y = 1.0;
    return y

@jit(nopython=True)
def fznd(d_formation, z_formation):
    d_spread = 3  # How far we will move from d_formation on each side
    steps_in_each_direction = 300 # How many steps d_spread will be divided into
    #In principle, the next 3 variables should be calculated dynamically, but for now I'm just fixing them for performance testing purposes
    step_length = 0.01  # = d_spread/steps_in_each_direction
    half_step_length = 0.005  # = 0.5*step_length
    total_steps = 601 # = 1 + (2*steps_in_each_direction)
    x = np.zeros(total_steps)
    i = d_formation - d_spread
    prev_normcdf = normcdf(i-half_step_length, d_formation, z_formation)
    for s in frange(0, steps_in_each_direction + 1, 1):
        new_normcdf = normcdf(i+half_step_length, d_formation, z_formation)
        x[s] = new_normcdf - prev_normcdf
        prev_normcdf = new_normcdf
        i += step_length  
    x[301:601] = x[0:300][::-1]
    #for s in frange(steps_in_each_direction + 1, total_steps, 1):
    #    x[s] = x[(total_steps - 1) - s]
    return x
