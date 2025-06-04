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

dt_r2_1D = np.array([ .009, .0595,  .2635,  2.012, 7.1435])
n_r2_1D  = np.array([2**20, 2**22,  2**24,  2**26,  2**28])
dt_r4_1D = np.array([.0095,  .058,   .264, 2.0005,  7.201])
n_r4_1D  = np.array([2**20, 2**22,  2**24,  2**26,  2**28])


dt_r2_2D = np.array([ .001,  .006,  .0565,   .290,  1.360])


n_r2_2D  = np.array([2**20, 2**22,  2**24,  2**26,  2**28])
dt_r4_2D = np.array([ .001,         .0605,         1.2625])
n_r4_2D  = np.array([2**20,         2**24,          2**28])
dt_r2_3D = np.array([.0025,         .0225,           .191])
n_r2_3D  = np.array([2**21,         2**24,          2**27])
dt_r4_3D = np.array([               .0225                ])
n_r4_3D  = np.array([               2**24                ])

dt_per_el_r2_1D = dt_r2_1D / n_r2_1D
dt_per_el_r4_1D = dt_r4_1D / n_r4_1D
dt_per_el_r2_2D = dt_r2_2D / n_r2_2D
dt_per_el_r4_2D = dt_r4_2D / n_r4_2D
dt_per_el_r2_3D = dt_r2_3D / n_r2_3D
dt_per_el_r4_3D = dt_r4_3D / n_r4_3D

fontsize = 15

fig, ax = plt.subplots(nrows=1, ncols=2, sharey=True, figsize=(8, 4.5))
ax[0].loglog(n_r2_1D, dt_per_el_r2_1D, 'o-', label="1D", color="b")
ax[0].loglog(n_r2_2D, dt_per_el_r2_2D, 'o-', label="2D", color="r")
ax[0].loglog(n_r2_3D, dt_per_el_r2_3D, 'o-', label="3D", color="g")
ax[1].loglog(n_r4_1D, dt_per_el_r4_1D, 'o-', label="1D", color="b")
ax[1].loglog(n_r4_2D, dt_per_el_r4_2D, 'o-', label="2D", color="r")
ax[1].loglog(n_r4_3D, dt_per_el_r4_3D, 'o-', label="3D", color="g")
ax[0].set_ylabel("Per-element permutation cost, "  + r'$s$', fontsize=fontsize)
ax[0].set_xlabel(r'$n=n_1\cdots n_d$', fontsize=fontsize)
ax[1].set_xlabel(r'$n=n_1\cdots n_d$', fontsize=fontsize)
ax[0].grid(), ax[1].grid()
ax[0].tick_params(axis='both', which='major', labelsize=fontsize)
ax[1].tick_params(axis='both', which='major', labelsize=fontsize)
ax[0].set_title("Radix-2", fontsize=fontsize)
ax[1].set_title("Radix-4", fontsize=fontsize)
plt.legend(fontsize=fontsize)
plt.savefig("per-element-permutation-cost.pdf")