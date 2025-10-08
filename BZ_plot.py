# -*- coding: utf-8 -*-
"""
Created on Fri Feb 11 13:29:55 2022

@author: Sven
"""

from abipy.abilab import abiopen
import matplotlib.pyplot as plt
from matplotlib import cm
import numpy as np
import pymatgen.electronic_structure.plotter as pl


def plot_BZ(bz_lattice, kpoints=None, ax=None, **kwargs):
    import pymatgen.electronic_structure.plotter as pl
    fig = None
    if ax is None:
        fig, ax = pl.plot_lattice_vectors(bz_lattice, ax=ax, color='red')
        pl.plot_wigner_seitz(bz_lattice, ax=ax)
    ax.set_xlim3d(-1, 1)
    ax.set_ylim3d(-1, 1)
    ax.set_zlim3d(-1, 1)

    if kpoints is not None:
        if kpoints.shape[0] > 4:
            # print('multiple')
            for k in kpoints:
                x, y, z = k[:3]
                ax.scatter(x, y, z, **kwargs)
        else:
            # print('single')
            x, y, z = kpoints[:3]
            ax.scatter(x, y, z, **kwargs)

    # ax.set_aspect('equal')
    ax.axis("off")

    return fig


# Path to the directory containing the GSR.nc file and Ncm.csv file
# path = "C:\\Users\\Sven\\OneDrive\\Studie\\Master\\Jaar 2\\Thesis\\
# Python scripts\\Data\\MoTe2\\10x10x5_0.90eV_30Ha_Emin"
# dirname = path[-path[::-1].index('\\'):]
path = "C:\\Users\\sbranchett\\work\\sven"
dirname = path[-path[::-1].index('\\'):]
# Name of the GSR.nc file
# abinit_file = path+"\\bulk_1o_GSR.nc"
abinit_file = "MoTe2_4x4x4_1o_GSR.nc"

# Color the k-points
give_color = True

# Size of circles representing the kpoints
s = 48

# Set to True if you want a top view or side view figure
top = True
side = True

with abiopen(abinit_file) as ncfile:
    abinit = ncfile.ebands

bz_lattice = abinit.kpoints.reciprocal_lattice
kpoints = abinit.kpoints.get_cart_coords()

# correcting kpoints from reciprocal lattice unit cell to first brillouin zone
reciprocal_lattice = abinit.kpoints.reciprocal_lattice.matrix
a = np.linalg.norm(reciprocal_lattice[0])

################################################
# The following part only works for hexagonal lattices.
# Please remove if other lattice is used
for ik, k in enumerate(kpoints):
    x, y, z = k
    if y <= np.sqrt(3)*x - a + 0.01 and y < 0:
        kpoints[ik] = k - reciprocal_lattice[0]
    if y >= - np.sqrt(3)*x + a - 0.01 and y >= 0:
        kpoints[ik] = k - reciprocal_lattice[1]
    if y < - np.sqrt(3)*x - a and y < 0:
        kpoints[ik] = k + reciprocal_lattice[1]
    if y > np.sqrt(3)*x + a and y > 0:
        kpoints[ik] = k + reciprocal_lattice[0]
    if y == 0 and y > np.sqrt(3)*x + a and x < 0:
        kpoints[ik] = k + reciprocal_lattice[0]
################################################


Ncm = np.genfromtxt(path+"\\Ncm.csv", delimiter=',') / abinit.nkpt**2
ksum = [sum(k) for k in Ncm]
norm = cm.colors.Normalize(vmin=min(ksum), vmax=max(ksum))
# norm = [np.amin(Ncm), np.amax(Ncm)]
viridis = plt.colormaps["viridis"]

ksum = [sum(k) for k in Ncm]
ksum = (ksum - min(ksum))
ksum = ksum / max(ksum)

# testing color map
# x = np.linspace(0,1,len(cmap))
# y = np.linspace(0,256,len(cmap))
# for k in range(len(cmap)):
#     plt.scatter(x[k], y[k], color=viridis(1.1))
# print(max(cmap))
# print(min(cmap))

if top:
    fig, ax = pl.plot_lattice_vectors(bz_lattice, color='red')
    for k in range(len(kpoints)):
        if give_color:
            color = viridis(ksum[k])
            plot_BZ(
                bz_lattice, kpoints[k], ax=ax, s=s, color=viridis(ksum[k])
            )
        else:
            pl.plot_BZ(bz_lattice, kpoints[k], ax=ax, s=s, color='b')
    pl.plot_wigner_seitz(bz_lattice, ax=ax)
    ax.view_init(elev=90, azim=-90)
    cbar = fig.colorbar(
        cm.ScalarMappable(norm=norm, cmap=viridis), ax=ax,
        shrink=0.5, label=r'$N_{CM} /; / /; {n_k}^2$'
    )
    cbar.set_label(r'$N_{CM} \; / \; {n_k}^2$', fontsize=12)
    cbar.ax.tick_params(labelsize=12)
    ax.text(0, 1.2, 1, dirname)
    fig.savefig(dirname+"_top.png", dpi=300)
    plt.show()

if side:
    fig, ax = pl.plot_lattice_vectors(bz_lattice, color='red')
    for k in range(len(kpoints)):
        if give_color:
            color = viridis(ksum[k])
            plot_BZ(bz_lattice, kpoints[k], ax=ax, s=s, color=viridis(ksum[k]))
        else:
            plot_BZ(bz_lattice, kpoints[k], ax=ax, s=s, color='b')
    pl.plot_wigner_seitz(bz_lattice, ax=ax)
    cbar = fig.colorbar(
        cm.ScalarMappable(norm=norm, cmap=viridis), ax=ax, shrink=0.5,
        label=r'$N_{CM} /; / /; {n_k}^2$'
    )
    cbar.set_label(r'$N_{CM} \; / \; {n_k}^2$', fontsize=12)
    ax.text(-0.63, 0.8, 1.1, dirname)
    ax.view_init(elev=43, azim=-75)
    fig.savefig(dirname+"_side.png", dpi=300)
    plt.show()
