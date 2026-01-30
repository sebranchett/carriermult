# -*- coding: utf-8 -*-
"""
Created on Wed Jan  5 19:20:18 2022

@author: Sven
"""

import numpy as np
import pandas as pd
import time
import os
from os import path
import sys
from datetime import datetime
from datetime import timedelta
from abipy.abilab import abiopen
from yambopy import YamboElectronsDB, YamboLatticeDB, YamboQPDB
import multiprocessing
from tqdm import tqdm
from joblib import Parallel, delayed
import argparse
import yaml

###############################################################################
# Function definitions
###############################################################################


def msg_title(msg):
    n = (75 - 2 - len(msg))//2
    log("#"*n + ' ' + msg + ' ' + "#"*n)


def msg_warning(msg):
    log()
    log("#"*75)
    log('# Warning:')
    log('# ' + msg)
    log("#"*75)
    log()


def msg_error(msg):
    log()
    log("#"*75)
    log('# Error:')
    log('# ' + msg)
    log("#"*75)
    log()


def log(info=''):
    print(info)
    with open(logfilename, 'a') as logfile:
        logfile.write(info)
        logfile.write('\n')


def print_progress(i, imax, start_time, filename):
    progress = i*100/(imax)
    filesize = path.getsize(filename)/1E6
    pred_filesize = filesize / (progress / 100)
    runtime = time.time() - start_time
    log('i: {}/{}\t ({}%)\nPredicted size: {:1.3f} MB'.format(
        i, imax, round(progress, 1), pred_filesize
    ))
    log('\nRuntime: {} \t Done in: {}'.format(
        timedelta(seconds=runtime),
        timedelta(seconds=(runtime)/(progress/100)-runtime)
    ))
    log('-'*75)


def check_moment(i, f, ii, ff, tol=1E-6, debug=False):
    net = abs(
        kpoints[int(i)]-kpoints[int(f)]+kpoints[int(ii)]-kpoints[int(ff)]
    )[:3]
    if debug:
        log(
            f"\ti: {kpoints[int(i)]} "
            f"\n-\tf: {kpoints[int(f)]} "
            f"\n+\tii: {kpoints[int(ii)]} "
            f"\n-\tff: {kpoints[int(ff)]} = {net}"
        )
    if (net <= np.array([tol, tol, tol])).all():
        return True
    else:
        return False


def check_IBZ(kpoint):
    if 0 <= kpoint[0] <= 0.5*reciprocal_lattice[0, 0]:
        if (
            0.5*reciprocal_lattice[0, 1] <= kpoint[1] <=
            0.5*reciprocal_lattice[1, 1]
        ):
            if 0 <= kpoint[2] <= 0.5*reciprocal_lattice[2, 2]:
                return True


def k2frac(kpoint):
    return np.dot(np.transpose(reciprocal_lattice_inv), kpoint[:3])


def frac2k(kpoint):
    return np.dot(np.transpose(reciprocal_lattice), kpoint[:3])


def red2e(bandindex):
    global usebands, degenerate
    if degenerate is True:
        return (2*bandindex + usebands[0]+1)
    if degenerate is False:
        return (bandindex + usebands[0]+1)


def e2red(bandindex):
    global usebands, degenerate
    if bandindex < usebands[0] + 1:
        raise ValueError(
            "Bandindex cannot be less than the lowest band used for "
            "red_energies"
        )
    if bandindex >= usebands[1]:
        raise ValueError(
            "Bandindex cannot be larger than the highest band used for "
            "red_energies"
        )
    if degenerate is True:
        return (bandindex - usebands[0] - 1)//2


def Umklapp(k):
    # I think this should work for all lattices, but I have not tested it for
    # anything other than hexagonal
    # Takes a kpoint (Cartesian coordinates) as input and converts that point
    # to a value within |0.5*reciprocal lattice vector|
    k_new = np.copy(k)
    k_new = k2frac(k_new[:3])

    for i, kx in enumerate(k_new):
        if abs(kx) > 0.5:
            k_new[i] = kx + 1 * np.sign(kx) * (np.sign(kx) * -1 *
                                               ((kx - 0.5) // (1) + 1))
        if abs(k_new[i] + 0.5) < 1E-6:
            k_new[i] = 0.5

    k_new = frac2k(k_new)
    return k_new


def find_CM_transitions(krow, energies, BG, Emin, Emax, save_CMfile=True,
                        etol=0.1):
    Ncm = np.zeros(energies.shape)
    # Spin character should still be added in this function!
    ki, kf, kii, kff = krow[:]
    if save_CMfile is True:
        transitions = np.zeros([0, 4])
    # transitions = np.zeros([0,8]) # save transitions as integer indices
    # Filter Ei (initial state at ki)
    Ei_vals = energies[ki]
    valid_Ei = np.where((Ei_vals >= Emin) & (Ei_vals <= Emax))[0]
    Eii_vals = energies[kii]
    valid_Eii = np.where((Eii_vals <= 0))[0]
    Ef_vals = energies[kf]
    # valid depends on Ei, so will be defined later
    Eff_vals = energies[kff]
    valid_Eff = np.where((Eff_vals > 0))[0]
    for Ei in valid_Ei[::-1]:
        Ei_energy = Ei_vals[Ei]
        valid_Ef = np.where((Ef_vals > Emin) & (Ei_energy - Ef_vals >= BG))[0]
        for Ef in valid_Ef[::-1]:
            Ef_energy = Ef_vals[Ef]
            if Ei_energy > 0 and Ef_energy <= 0:
                # Only consider full VB or full CB transitions. If electron
                # goes from CB to VB then there is no CM).
                continue
            Eif_energy = Ei_energy - Ef_energy
            for Eii in valid_Eii:
                Eifii_energy = Eif_energy + Eii_vals[Eii]
                for Eff in valid_Eff[::-1]:
                    dE = Eifii_energy - Eff_vals[Eff]
                    if abs(dE) < etol:
                        # For CB the initial CM state is an electron
                        if Ei_energy > 0:
                            new_transition = [[ki, Ei, kf, Ef]]
                        # For a VB transition, the CM carrier is a hole and
                        # therefore we change the order of initial & final
                        # state
                        else:  # if Ei_energy <= 0:
                            new_transition = [[kf, Ef, ki, Ei]]
                        Ncm[new_transition[0][0], new_transition[0][1]] += 1
                        if save_CMfile is True:
                            transitions = np.append(
                                transitions, new_transition, axis=0
                            )
                        # Need to add a part where it saves truth values for
                        # all energies
    if save_CMfile is True:
        if len(transitions) == 0:
            transitions = [[None]]
        return Ncm, transitions
    else:
        return Ncm


def calculate_CM_transitions(
    i, kpoints, energies, BG, Emin, Emax, save_kfile=False, kfilename=None,
    save_CMfile=False, cmfilename=None, chunksize=10000, ktol=1E-6, etol=0.01
):
    global count, start_time
    k_range = (kpoints.shape[0])
    Ncm = np.zeros(energies.shape)
    ktransitions = np.zeros([0, 4], dtype=int)
    cmtransitions = np.zeros([0, 4], dtype=int)
    for f in range(k_range):
        for ii in range(k_range):
            kff = kpoints[i, :] + kpoints[ii, :] - kpoints[f, :]
            kff = Umklapp(kff)
            new_length = (kff[0]**2 + kff[1]**2 + kff[2]**2)**0.5
            kff = np.append(kff, new_length)
            under = kff[-1] - ktol
            above = kff[-1] + ktol
            a = np.where(
                np.logical_and(
                    kpoints[:, -1] >= under,
                    kpoints[:, -1] <= above
                )
            )
            if a[0].shape[0] == 0:
                print(
                    f'{k2frac(kpoints[i, :])} + {k2frac(kpoints[ii, :])} - '
                    f'{k2frac(kpoints[f, :])} = {k2frac(kff)}'
                )
                # print(f'Mismatching k-point: {k2frac(kff)}\n')
            for ff in a[0]:
                # this may be improved with a np.find function. Find index of
                # values between range
                # if kpoints[ff, -1] - kff[-1] < ktol:
                # if not check_moment(i,f,ii,ff,ktol):
                if not (kpoints[ff, :3] - kff[:3] < 1E-6).all():
                    if ff == a[0][-1]:
                        print(
                            f'{k2frac(kpoints[i, :])} + '
                            f'{k2frac(kpoints[ii, :])} - '
                            f'{k2frac(kpoints[f, :])} = '
                            f'{k2frac(kff)}'
                        )
                        print('Warning: No matching k value found')
                    continue
                else:
                    ktransition = np.array([[i, f, ii, ff]])

                    if save_kfile:
                        ktransitions = np.append(
                            ktransitions, ktransition, axis=0
                        )
                        if len(ktransition) > chunksize:
                            with open(kfilename, 'a') as karray_csv:
                                np.savetxt(
                                    karray_csv, ktransitions, delimiter=',',
                                    fmt='%i'
                                )
                            ktransitions = np.zeros([0, 4], dtype=int)

                    if save_CMfile:
                        Ncm_in, cm = find_CM_transitions(
                            ktransition[0], energies, BG, Emin, Emax,
                            save_CMfile=save_CMfile, etol=etol
                        )
                        if not (np.array(cm) is None).any():
                            cmtransitions = np.append(
                                cmtransitions, cm, axis=0
                            )
                        if len(cm) > chunksize or len(ktransition) > chunksize:
                            with open(cmfilename, 'a') as cm_csv:
                                np.savetxt(
                                    cm_csv, cmtransitions, delimiter=',',
                                    fmt='%i'
                                )
                            cmtransitions = np.zeros([0, 4], dtype=int)

                    if not save_CMfile:
                        Ncm_in = find_CM_transitions(
                            ktransition[0], energies, BG, Emin, Emax,
                            save_CMfile=save_CMfile, etol=etol
                        )

                    Ncm += Ncm_in
                    break

    if save_CMfile:
        with open(cmfilename, 'a') as cm_csv:
            np.savetxt(cm_csv, cmtransitions, delimiter=',', fmt='%i')
        cmtransitions = np.zeros([0, 4], dtype=int)

    if save_kfile:
        with open(kfilename, 'a') as karray_csv:
            np.savetxt(karray_csv, ktransitions, delimiter=',', fmt='%i')
        ktransitions = np.zeros([0, 4], dtype=int)

    return Ncm


def read_CM_states(CMfilename, energies):
    Ncm = np.zeros(energies.shape)
    CMs = pd.read_csv(CMfilename, names=['ki', 'ni', 'kf', 'nf'])
    for order in range(1):
        if order == 0:
            for index, row in CMs.iterrows():
                ki, ni = int(row['ki']), int(row['ni'])
                Ncm[ki, ni] += 1


def plot_BZ(bz_lattice, kpoints=None, ax=None, **kwargs):
    import pymatgen.electronic_structure.plotter as pl
    fig = None
    if ax is None:
        fig, ax = pl.plot_lattice_vectors(bz_lattice, ax=ax)
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


def find_CMcount(
    energies, Ncm, Emax=4, dE=0.01, degenerate=True, averaging=False
):
    # returns array with
    # col0: photon energy that will produce electron-hole pair
    # col1: # primary electron CM transitions
    # col2: # primary hole CM transitions
    if degenerate:
        factor = 2  # This is needed because Ncm was accounted for degeneracy
        #             by multiplying with 16
        #             But if we excite at a given energy, then the electron &
        #             hole state has 2 bands it can occupy
        #             The number of transitions that are possible from these
        #             degenerate bands is two times smaller
        #             than if we kept viewing it as a single band. In other
        #             word, we counted all transitions as if
        #             all degenerate bands were a single band, but if we now
        #             consider them two seperate bands again,
        #             naturally we should divide the number of transitions by
        #             2.
    CMmatrix = np.zeros([0, 3])
    for k, Ek in enumerate(tqdm(energies)):
        for Eh in range(len(Ek)):  # Eh: bandindex of hole
            if energies[k, Eh] > 0:
                break
            for Ee in range(len(Ek)):  # Ee: bandindex of electron
                if energies[k, Ee] <= 0:
                    continue
                Ephoton = energies[k, Ee] - energies[k, Eh]
                if Ephoton > Emax:
                    continue
                Ncme, Ncmh = Ncm[k, Ee]/factor, Ncm[k, Eh]/factor
                # number of CM first order transitions for electron and hole

                CMmatrix = np.append(
                    CMmatrix, np.array([[Ephoton, Ncme, Ncmh]]), axis=0
                )
    CMmatrix = CMmatrix[np.argsort(CMmatrix[:, 0])]

    if averaging:
        df = pd.DataFrame(CMmatrix, columns=['E', 'Ne', 'Nh'])
        E = np.arange(0, 4, dE)
        Ncme = []
        Ncmh = []

        for ix, Ex in enumerate(E):
            if ix == 0:
                Ncme.append(0)
                Ncmh.append(0)
                continue
            erange = (Ex-dE/2, Ex+dE/2)
            Ncme_mean = df[df['E'].between(*erange)]['Ne'].mean()
            Ncmh_mean = df[df['E'].between(*erange)]['Nh'].mean()
            if not Ncme_mean >= 0:
                Ncme_mean = 0
            if not Ncmh_mean >= 0:
                Ncmh_mean = 0
            Ncme.append(Ncme_mean)
            Ncmh.append(Ncmh_mean)
        CMmatrix = np.zeros([len(E), 3])
        CMmatrix[:, 0], CMmatrix[:, 1], CMmatrix[:, 2] = E, Ncme, Ncmh

    return CMmatrix


def N2QY(N, rate):
    # Convert number of CM states to quantum yield with paramater rate. rate
    # stands for cooling rate to carrier multiplication ratio
    # rate = Rcool/Rcm
    return N / (N + rate)


def calculate_quantum_yield(CMmatrix, Re, Rh):
    QY = np.zeros([len(CMmatrix), 2])
    QY[:, 0] = CMmatrix[:, 0]
    for i in range(len(QY)):
        QY[i, 1] = 1 + N2QY(CMmatrix[i, 1], Re) + N2QY(CMmatrix[i, 2], Rh)
    return QY


def load_and_organize_CMcsv(CMfilename):
    # Loads the CM transitions file as a pandas dataframe, adjusts initial
    # states in valance band (initial state = final state)
    # and sorts dataframe by initial k and E
    CM_transitions = pd.read_csv(
        CMfilename, names=['ki', 'Ei', 'kf', 'Ef']
    ).sort_values(by=['Ei'])
    for index, row in CM_transitions.iterrows():
        if row['Ei'] > 0:
            break
        Ei, Ef = row['Ei'], row['Ef']
        ki, kf = row['ki'], row['kf']
        CM_transitions.iloc[index, 0] = int(kf)
        CM_transitions.iloc[index, 1] = Ef
        CM_transitions.iloc[index, 2] = int(ki)
        CM_transitions.iloc[index, 3] = Ei
    return CM_transitions.sort_values(by=['ki', 'Ei', 'kf', 'Ef'])


def load_abinit_nc_file(abinit_file):
    log()
    msg_title('INITIALIZING INPUT DATA')
    log()
    log("Loading abinit file from: \n"+abinit_file)
    log()
    try:
        with abiopen(abinit_file) as ncfile:
            abinit = ncfile.ebands
    except FileNotFoundError:
        msg_error(
            'Abinit *GSR.nc file cannot be loaded, '
            'check given path for Abinit file'
        )
        sys.exit()

    log("File loaded")
    info = [abinit.nkpt, abinit.nband]
    log("k-points: {} \nBands: {}".format(*info))
    log()

#   ############################################

#   ############### gather kpoints #############

    log(
        "Gathering k-points in Cartesian coordinates & "
        "calculating vector magnitudes.."
    )
    reciprocal_lattice = abinit.kpoints._reciprocal_lattice._matrix
    reciprocal_lattice_inv = abinit.kpoints._reciprocal_lattice.inv_matrix
    nkpt = abinit.nkpt
    kpoints = np.zeros([nkpt, 4])
    kpoints[:, :3] = abinit.kpoints.get_cart_coords()
    # The magnitude of the k-vectors are calculated in order to speed up the
    # find_ktransitions function
    for k in range(kpoints.shape[0]):
        kpoints[k, 3] = np.sqrt(np.dot(kpoints[k, :3], kpoints[k, :3]))

    log("Done\n")

#   ############################################

#   ############# gather eigenvalues ###########
    log("Gathering energies..")
    energies = (np.copy(abinit._eigens[0]) - abinit.get_e0('fermie')).to('eV')
    # sets VBM as 0 and converts to energies to eV
    # Since the BG is underestimated in DFT, a scissor operator is used to
    # shift the CB so that it matches with experimental data
    CALC_BG = abinit.fundamental_gaps[0].energy
    scissor = TrueBG - CALC_BG
    log('Scissor operation:')
    log(f'Calculated bandgap: {CALC_BG} eV')
    log(f'Experimental bandgap: {TrueBG} eV')
    log(f'Shifting CB with: {scissor:1.6f} eV')
    log('Please check if the calculated bandgap is indeed in eV.')
    log()

    cbm = abinit.lumos[0][2]
    energies[:, cbm:] = energies[:, cbm:] + scissor

    return reciprocal_lattice, reciprocal_lattice_inv, \
        kpoints, energies


def load_yambo_nc_file(
        yambo_dir, yambo_file='ns.db1', gw_dir='', gw_file='ndb.QP'):
    log()
    msg_title('INITIALIZING INPUT DATA')
    log()
    log("Loading Yambo file from: \n" + yambo_dir + "/" + yambo_file)
    log()
    try:
        yambo = YamboElectronsDB.from_db_file(yambo_dir, yambo_file)
        lattice = YamboLatticeDB.from_db_file(yambo_dir + "/" + yambo_file)
        nbands = yambo.nbands
    except FileNotFoundError:
        msg_error(
            'Yambo ns.db1 file cannot be loaded, '
            'check given path and filename'
        )
        sys.exit()
    if gw_dir != '':
        try:
            ydb = YamboQPDB.from_db(folder=gw_dir, filename=gw_file)
            nbands = ydb.nbands
        except FileNotFoundError:
            msg_error(
                'Yambo GW ndb.QP file cannot be loaded, '
                'check given path and filename'
            )
            sys.exit()

    log("File loaded")
    nkpt = lattice.nkpoints
    info = [nkpt, nbands]
    log("k-points: {} \nBands: {}".format(*info))
    log()

    #   ############################################

    #   ############### gather kpoints #############

    log(
        "Gathering k-points in Cartesian coordinates & "
        "calculating vector magnitudes.."
        )

    # The Umklapp function needs the reciprocal lattice and its inverse
    global reciprocal_lattice, reciprocal_lattice_inv
    reciprocal_lattice = lattice.rlat
    reciprocal_lattice_inv = np.linalg.inv(reciprocal_lattice)
    kpoints = np.zeros([nkpt, 4])
    kpoints[:, :3] = lattice.car_kpoints

    # The Yambo expansion of kpoints from the IBZ leads to
    # not all kpoints in the first BZ
    for i in range(kpoints.shape[0]):
        kpoints[i, :3] = Umklapp(kpoints[i, :3])
        # The magnitude of the k-vectors are calculated in order to speed
        # up the find_ktransitions function
        kpoints[i, 3] = np.sqrt(np.dot(kpoints[i, :3], kpoints[i, :3]))

    log("Done\n")

    #   ############################################

    #   ############# gather eigenvalues ###########
    log("Gathering energies..")

    # Expand the eigenvalues from the IBZ to the full BZ
    yambo.expandEigenvalues()
    cbm = yambo.nbandsv

    if not gw_dir:
        # Convert to np.float64 for faster results
        energies = (yambo.eigenvalues[0, :, :]).astype(np.float64, copy=False)
    else:
        cbm = cbm-ydb.min_band+1  # Not all bands are calculated in GW
        # energies = ydb.expand_eigenvalues(
        #     lattice=lattice, data=ydb.eigenvalues_dft)
        energies = ydb.expand_eigenvalues(
            lattice=lattice, data=ydb.eigenvalues_qp)

    max_valence_energy = np.max(energies[:, :cbm])
    min_conduction_energy = np.min(energies[:, cbm:])
    # sets VBM as 0
    # Since the BG is underestimated in DFT, a scissor operator is used to
    # shift the CB so that it matches with experimental data
    energies = energies - max_valence_energy
    CALC_BG = min_conduction_energy - max_valence_energy
    scissor = TrueBG - CALC_BG
    log('Scissor operation:')
    log(f'Calculated bandgap: {CALC_BG} eV')
    log(f'Experimental bandgap: {TrueBG} eV')
    log(f'Shifting CB with: {scissor:1.6f} eV')
    log('Please check if the calculated bandgap is indeed in eV.')
    log()

    energies[:, cbm:] = energies[:, cbm:] + scissor

    return reciprocal_lattice, reciprocal_lattice_inv, \
        kpoints, energies


###############################################################################
# Data Initialization
###############################################################################
if __name__ == "__main__":      # This is needed if you want to import
    #                         functions from this script to another python file
    # Parse command-line arguments to optionally override input files
    parser = argparse.ArgumentParser(
        description=(
            'Calculate carrier multiplication transitions from Abinit '
            'or Quantum ESPRESSO/Yambo data.'
        )
    )
    parser.add_argument(
        '-f',
        '--file',
        help='Path to yaml file containing input parameters',
        metavar='INPUT_FILE',
        default='cm_input.yaml',
    )
    args = parser.parse_args()
    # Read input from yaml file
    with open(args.file, 'r') as infile:
        input_data = yaml.safe_load(infile)
    abinit_file = input_data.get('abinit_file', '')
    yambo_dir = input_data.get('yambo_dir', '')
    yambo_file = input_data.get('yambo_file', 'ns.db1')
    gw_dir = input_data.get('gw_dir', '')
    gw_file = input_data.get('gw_file', 'ndb.QP')
    TrueBG = input_data.get('TrueBG')  # eV
    Emax = input_data.get('Emax')  # eV
    etol = input_data.get('etol')  # eV
    run_CM_calculations = input_data.get('run_CM_calculations', True)
    save_kfile = input_data.get('save_kfile', False)
    kfilename = input_data.get('kfilename', 'k_transitions.csv')
    save_CMfile = input_data.get('save_CMfile', False)
    CMfilename = input_data.get('CMfilename', 'CM_transitions.csv')
    remove_degenerate = input_data.get('remove_degenerate', True)
    logfilename = input_data.get('logfilename', 'cm_logfile.log')

    if path.exists(logfilename):
        os.remove(logfilename)
    log('\n\n')
    log('File made on {}'.format(datetime.now().strftime("%Y-%m-%d %H:%M:%S")))
    msg_title('INPUT')
    log()
    log('run_CM_calculations:\t\t{}'.format(run_CM_calculations))
    log()
    if save_kfile:
        log('kfilename:\t\t{}'.format(kfilename))
    log()
    if save_CMfile:
        log('\nCMfilename:\t\t{}'.format(CMfilename))
    log(f'save_kfile\t\t\t{save_kfile}\nsave_CMfile\t\t\t{save_CMfile}')
    log()
    Emin = TrueBG - Emax
    log(
        f"Emin:\t\t\t\t{Emin} eV\n"
        f"Emax:\t\t\t\t{Emax} eV\n"
        f"Experimental Eg:\t{TrueBG} eV\n"
        f"Error tolerance\t\t{etol} eV"
    )
    # log('Limiting energy calculations to {}% of the {} data'.
    # format(datasize, kfilename))
    log(f'remove_degenerate\t\t{remove_degenerate}')
    log()
    log(f'Calculating with {multiprocessing.cpu_count()} cores')
    msg_title('END INPUT')

    # ############## Load dft file ###############
    # if abinit file has been specified:
    if abinit_file:
        # check abinit_file exists
        if not path.exists(abinit_file):
            msg_error(
                f'Abinit file {abinit_file} not found. Please check path.'
            )
            sys.exit()
        reciprocal_lattice, reciprocal_lattice_inv, \
            kpoints, energies = load_abinit_nc_file(abinit_file)
    else:
        # check yambo_dir/yambo_file exists
        yambo_file_path = path.join(yambo_dir, yambo_file)
        if not path.exists(yambo_file_path):
            msg_error(
                f'Yambo file {yambo_file_path} not found. Please check path.'
            )
            sys.exit()
        # if gw_dir and gw_file have been specified, check they exist
        if gw_dir and gw_file:
            if not path.exists(path.join(gw_dir, gw_file)):
                msg_error(
                    f'Yambo GW file {path.join(gw_dir, gw_file)} not found. '
                    'Please check path.'
                )
                sys.exit()
        reciprocal_lattice, reciprocal_lattice_inv, \
            kpoints, energies = load_yambo_nc_file(
                yambo_dir, yambo_file, gw_dir, gw_file)
    #####
    # The following part limits the energy array to the values Emin and Emax
    nmin = 0
    nmax = 0
    _emax = 0
    for n in range(energies.shape[1]):
        if max(energies[:, n]) > _emax:
            _emax = max(energies[:, n])
        if max(energies[:, n]) >= Emin and nmin == 0:
            nmin = n
        if min(energies[:, n]) > Emax:
            nmax = n
            break
    if nmax == 0:
        nmax = energies.shape[1] - 1
        Emax = _emax
        msg_warning(
            "Maximum energy is found to be {:1.4f} eV, which is lower than"
            " Emax.\nEmax will be adjusted to this value."
            .format(Emax)
        )
    usebands = (nmin, nmax)
    red_energies = np.copy(energies[:, usebands[0]:usebands[1]])
    log(
        "Limiting energy matrix: lowest band: {}. "
        "highest band: {}.".format(usebands[0], usebands[1])
    )
    #####
    # If remove_degenerate is true, this part will remove degenerate bands
    # from the calculations
    degenerate = False
    if remove_degenerate:
        log("Checking if all energies are degenerate")
        degenerate = True
        for i, Ek in enumerate(red_energies):
            for j, Ekn in enumerate(Ek):
                # if j is even, save energy, if odd compare to previous j
                if j % 2 == 0:
                    prev_energy = Ekn
                    continue
                else:
                    if Ekn - prev_energy > 0.0001:
                        degenerate = False
                        log("Some energies are found to be non-degenerate")
                        break
            if not degenerate:
                break

        if degenerate:
            log(
                "All energies are found to be degenerate. "
                "Removing half of the energies for faster computation"
            )  # maybe need to add these again at some point for QY calculation
            red_energies = np.delete(
                red_energies, list(np.arange(0, red_energies.shape[1], 2)), 1
            )
    np.savetxt('reduced_energies.csv', red_energies)
    log("Done. Reduced energy array saved in \'reduced_energies.csv\'")
    log(
        "Note: to convert reduced band indices back to real bandindices, "
        "use the red2e(bandindex) function."
    )
    log()
    msg_title('END INITIALIZING INPUT DATA')
    log()

    ############################################
    # Define precision for output files
    float_precision = f'%1.{len(str(etol))-1}f'

    ###########################################################################
    # Run the actual calculations
    ###########################################################################

    if run_CM_calculations is True:
        msg_title('STARTING CALCULATIONS')
        log()
        if save_kfile:
            if path.exists(kfilename):
                msg_error(
                    '{} already exists, please rename the output file before '
                    'running again'.format(kfilename)
                )
                sys.exit()

        if save_kfile:
            if path.exists(CMfilename):
                msg_error(
                    '{} already exists, please rename the output file before '
                    'running again'.format(CMfilename)
                )
                sys.exit()

        Ncm = np.zeros(red_energies.shape)
        kparralel = [i for i in range(len(kpoints))]
        num_cores = multiprocessing.cpu_count()
        chunksize = 10000
        start_time = time.time()

        if True:  # Run function in parallel mode
            data = Parallel(n_jobs=num_cores)(
                delayed(calculate_CM_transitions)(
                    i, kpoints, red_energies, TrueBG,
                    Emin=Emin, Emax=Emax, save_kfile=save_kfile,
                    kfilename=kfilename, save_CMfile=save_CMfile,
                    cmfilename=CMfilename, etol=etol
                ) for i in tqdm(kparralel))
            Ncm += sum(data)
        else:  # Run function in serial mode
            Ncm += calculate_CM_transitions(
                0, kpoints, red_energies, TrueBG,
                Emin=Emin, Emax=Emax, save_kfile=save_kfile,
                kfilename=kfilename, save_CMfile=save_CMfile,
                cmfilename=CMfilename, etol=etol)
            # msg_error('\'run_combined\' can only be run parallel')
            # sys.exit()
        if degenerate:
            Ncm *= 16   # Each band is 2fold degenerate. We have 4 energy
            #             levels that we compare
            #             when finding CM transitions. 2*2*2*2 = 16
        end_time = time.time()
        log('CM transition calculations done!')
        runtime = end_time-start_time
        log('Walltime: {} ()'.format(timedelta(seconds=runtime)))
        if save_kfile:
            log("k-point array data saved in {}".format(kfilename))
        if save_CMfile:
            log("CM transitions data saved in {}".format(CMfilename))
        with open('Ncm.csv', 'w') as ncm:
            np.savetxt(ncm, Ncm, delimiter=',', fmt='%i')
        log(
            "Number of transitions per state data saved in {}"
            .format('Ncm.csv')
        )
        log('#'*75)
        log()
