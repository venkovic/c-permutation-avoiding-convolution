"""
 *Copyright (c) 2025 Nicolas Venkovic
 * 
 *This file is part of c-permutation-avoiding-convolution.
 * 
 *This file is licensed under the MIT License.
 *For the full license text, see the LICENSE file in the root directory of this project.
"""

import numpy as np
import matplotlib.pyplot as plt

plt.rc('text.latex', preamble=r'\usepackage{amssymb} \usepackage{amsmath}')
plt.rc('text', usetex=True)
plt.rc('font', family='serif')

dt_r2_1D   = np.array([.0170, .0780,  .3320, 1.5445,   6.205])
dt_r2_1D_2 = np.array([.0140, .0650,  .5065, 2.5100, 11.1795])
n_r2_1D    = np.array([2**20, 2**22,  2**24,  2**26,   2**28])
dt_r4_1D   = np.array([.0135, .0535,   .249, 1.0815,  4.5890])
dt_r4_1D_2 = np.array([ .011, .0550,   .330,  1.563,  6.9505])
n_r4_1D    = np.array([2**20, 2**22,  2**24,  2**26,   2**28])
dt_r2_2D   = np.array([ .029, .2000, 1.0385,  4.591, 20.0615])
dt_r2_2D_2 = np.array([.0355, .1420, 1.1605,  5.829, 23.7825])
n_r2_2D    = np.array([2**20, 2**22,  2**24,  2**26,   2**28])
dt_r4_2D   = np.array([.0280,        1.0295,         19.8865])
dt_r4_2D_2 = np.array([.0340,        1.1345,         23.5950])
n_r4_2D    = np.array([2**20,         2**24,           2**28])
dt_r2_3D   = np.array([.0054,        0.7220,           6.528])
dt_r2_3D_2 = np.array([.0795,        0.8175,         10.0875])
n_r2_3D    = np.array([2**21,         2**24,           2**27])
dt_r4_3D   = np.array([              0.7035,                ])
dt_r4_3D_2 = np.array([              0.7785                 ])
n_r4_3D    = np.array([               2**24                 ])

dt_per_el_r2_1D = dt_r2_1D / n_r2_1D
dt_per_el_r4_1D = dt_r4_1D / n_r4_1D
dt_per_el_r2_2D = dt_r2_2D / n_r2_2D
dt_per_el_r4_2D = dt_r4_2D / n_r4_2D
dt_per_el_r2_3D = dt_r2_3D / n_r2_3D
dt_per_el_r4_3D = dt_r4_3D / n_r4_3D

dt_per_el_r2_1D_2 = dt_r2_1D_2 / n_r2_1D
dt_per_el_r4_1D_2 = dt_r4_1D_2 / n_r4_1D
dt_per_el_r2_2D_2 = dt_r2_2D_2 / n_r2_2D
dt_per_el_r4_2D_2 = dt_r4_2D_2 / n_r4_2D
dt_per_el_r2_3D_2 = dt_r2_3D_2 / n_r2_3D
dt_per_el_r4_3D_2 = dt_r4_3D_2 / n_r4_3D

fontsize = 15

fig, ax = plt.subplots(nrows=1, ncols=2, sharey=True, figsize=(8, 4.5))
ax[0].loglog(n_r2_1D, dt_per_el_r2_1D, 'o-', label="1D", color="b")
ax[0].loglog(n_r2_2D, dt_per_el_r2_2D, 'o-', label="2D", color="r")
ax[0].loglog(n_r2_3D, dt_per_el_r2_3D, 'o-', label="3D", color="g")
ax[1].loglog(n_r4_1D, dt_per_el_r4_1D, 'o-', label="1D", color="b")
ax[1].loglog(n_r4_2D, dt_per_el_r4_2D, 'o-', label="2D", color="r")
ax[1].loglog(n_r4_3D, dt_per_el_r4_3D, 'o-', label="3D", color="g")
ax[0].loglog(n_r2_1D, dt_per_el_r2_1D_2, 'o--', color="b")
ax[0].loglog(n_r2_2D, dt_per_el_r2_2D_2, 'o--', color="r")
ax[0].loglog(n_r2_3D, dt_per_el_r2_3D_2, 'o--', color="g")
ax[1].loglog(n_r4_1D, dt_per_el_r4_1D_2, 'o--', color="b")
ax[1].loglog(n_r4_2D, dt_per_el_r4_2D_2, 'o--', color="r")
ax[1].loglog(n_r4_3D, dt_per_el_r4_3D_2, 'o--', color="g")
ax[0].set_ylabel("Average per-element butterfly cost, "  + r'$s$', fontsize=fontsize)
ax[0].set_xlabel(r'$n=n_1\cdots n_d$', fontsize=fontsize)
ax[1].set_xlabel(r'$n=n_1\cdots n_d$', fontsize=fontsize)
ax[0].grid(), ax[1].grid()
ax[0].tick_params(axis='both', which='major', labelsize=fontsize)
ax[1].tick_params(axis='both', which='major', labelsize=fontsize)
ax[0].set_title("Radix-2", fontsize=fontsize)
ax[1].set_title("Radix-4", fontsize=fontsize)
plt.legend(fontsize=fontsize)
plt.savefig("per-element-butterfly-cost.pdf")