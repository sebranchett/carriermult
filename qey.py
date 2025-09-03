import numpy as np
from yambopy import YamboElectronsDB, YamboLatticeDB
import sys

logfilename = 'log.txt'
yambo_dir = "./MoTe2/SAVE"
yambo_file = "ns.db1"
TrueBG = 0.9


def log(info=''):
    print(info)
    with open(logfilename, 'a') as logfile:
        logfile.write(info)
        logfile.write('\n')


def msg_title(msg):
    n = (75 - 2 - len(msg))//2
    log("#"*n + ' ' + msg + ' ' + "#"*n)


def msg_error(msg):
    log()
    log("#"*75)
    log('# Error:')
    log('# ' + msg)
    log("#"*75)
    log()


def load_yambo_nc_file(yambo_dir, yambo_file):
    log()
    msg_title('INITIALIZING INPUT DATA')
    log()
    log("Loading Yambo file from: \n" + yambo_dir + "/" + yambo_file)
    log()
    try:
        yambo = YamboElectronsDB.from_db_file(yambo_dir, yambo_file)
        lattice = YamboLatticeDB.from_db_file(yambo_dir + "/" + yambo_file)
    except FileNotFoundError:
        msg_error(
            'Yambo ns.db1 file cannot be loaded, '
            'check given path for Yambo file'
        )
        sys.exit()

    log("File loaded")
    nkpt = lattice.nkpoints
    info = [nkpt, yambo.nbands]
    log("k-points: {} \nBands: {}".format(*info))
    log()

    #   ############################################

    #   ############### gather kpoints #############

    log(
        "Gathering k-points in Cartesian coordinates & "
        "calculating vector magnitudes.."
        )

    reciprocal_lattice = lattice.rlat
    reciprocal_lattice_inv = np.linalg.inv(reciprocal_lattice)

    kpoints = np.zeros([nkpt, 4])
    kpoints[:, :3] = lattice.car_kpoints
    # The magnitude of the k-vectors are calculated in order to speed up the
    # find_ktransitions function
    for k in range(kpoints.shape[0]):
        kpoints[k, 3] = np.sqrt(np.dot(kpoints[k, :3], kpoints[k, :3]))

    log("Done\n")

    #   ############################################

    #   ############# gather eigenvalues ###########
    log("Gathering energies..")

    # Expand the eigenvalues from the IBZ to the full BZ
    yambo.expandEigenvalues()

    cbm = yambo.nbandsv
    max_valence_energy = np.max(yambo.eigenvalues[0, :, :cbm])
    min_conduction_energy = np.min(yambo.eigenvalues[0, :, cbm:])
    energies = yambo.eigenvalues[0, :, :] - max_valence_energy
    # sets VBM as 0
    # Since the BG is underestimated in DFT, a scissor operator is used to
    # shift the CB so that it matches with experimental data
    DFT_BG = min_conduction_energy - max_valence_energy
    scissor = TrueBG - DFT_BG
    log('Scissor operation:')
    log(f'Calculated bandgap: {DFT_BG} eV')
    log(f'Experimental bandgap: {TrueBG} eV')
    log(f'Shifting CB with: {scissor:1.6f} eV')
    log('Please check if the calculated bandgap is indeed in eV.')
    log()

    energies[:, cbm:] = energies[:, cbm:] + scissor

    return reciprocal_lattice, reciprocal_lattice_inv, \
        kpoints, energies


load_yambo_nc_file(yambo_dir, yambo_file)
