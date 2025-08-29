import numpy as np
import yambopy as ymb

yambo_dir = "./MoTe2/SAVE"
yambo_file = "ns.db1"
yambo = ymb.YamboElectronsDB.from_db_file(yambo_dir, yambo_file)
lattice = ymb.YamboLatticeDB.from_db_file(yambo_dir + "/" + yambo_file)
yambo_eigs_ibz = yambo.eigenvalues_ibz

# print out the dimensions of yambo_eigs_ibz
print("Shape of IBZ eigenvalues:", yambo_eigs_ibz.shape)

yambo.expandEigenvalues()
yambo_eigs = yambo.eigenvalues
print("Shape of eigenvalues:", yambo_eigs.shape)
print("Number of occupied bands:", yambo.nbandsv)
print("Value at (0, 0, 91):", yambo_eigs[0][0][91])
print("Value at (0, 0, 92):", yambo_eigs[0][0][92])

print("Summary information about yambo_eigs:")
print("Shape:", yambo_eigs.shape)
print("Max:", np.max(yambo_eigs))
print("Min:", np.min(yambo_eigs))
print("Energy gaps:", yambo.energy_gaps())

max_valence = np.max(yambo_eigs[:, :, :92])
min_conduction = np.min(yambo_eigs[:, :, 92:])
print("Max valence band energy:", max_valence)
print("Min conduction band energy:", min_conduction)
print("Calculated band gap:", min_conduction - max_valence)

print("Reciprocal lattice:")
print(lattice.rlat)

nkpt = lattice.nkpoints
kpoints = np.zeros([nkpt, 4])
kpoints[:, :3] = lattice.car_kpoints
# The magnitude of the k-vectors are calculated in order to speed up the
# find_ktransitions function
for k in range(kpoints.shape[0]):
    kpoints[k, 3] = np.sqrt(np.dot(kpoints[k, :3], kpoints[k, :3]))

print("K points:")
print(kpoints)
