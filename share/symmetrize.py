# import ase.io
# import sys
import numpy as np
import spglib

def prep(at, symprec=1.0e-5):
    dataset = spglib.get_symmetry_dataset(at, symprec=symprec)
    print "symmetry.prep got symmetry group number",dataset["number"],", international (Hermann-Mauguin)",dataset["international"],", Hall",dataset["hall"]
    rotations = dataset['rotations'].copy()
    translations = dataset['translations'].copy()
    symm_map=[]
    scaled_pos = at.get_scaled_positions()
    for (r, t) in zip(rotations, translations):
        # print "op"
        this_op_map = [-1] * len(at)
        for i_at in range(len(at)):
            new_p = np.dot(r, scaled_pos[i_at,:]) + t
            dp = scaled_pos - new_p
            dp -= np.round(dp)
            i_at_map = np.argmin(np.linalg.norm(dp,  axis=1))
            this_op_map[i_at] = i_at_map
        symm_map.append(this_op_map)
    return (rotations, translations, symm_map)

# lattice vectors expected as row vectors (same as ASE get_cell() convention), inv_lattice is its matrix inverse (get_reciprocal_cell().T)
def forces(lattice, inv_lattice, forces, rot, trans, symm_map):
    scaled_symmetrized_forces_T = np.zeros(forces.T.shape)

    scaled_forces_T = np.dot(inv_lattice.T,forces.T)
    for (r, t, this_op_map) in zip(rot, trans, symm_map):
        # print "op "
        # print r, t, this_op_map
        transformed_forces_T = np.dot(r, scaled_forces_T)
        scaled_symmetrized_forces_T[:,this_op_map[:]] += transformed_forces_T[:,:]
    scaled_symmetrized_forces_T /= len(rot)

    symmetrized_forces = np.dot(lattice.T, scaled_symmetrized_forces_T).T

    return symmetrized_forces

# lattice vectors expected as row vectors (same as ASE get_cell() convention), inv_lattice is its matrix inverse (get_reciprocal_cell().T)
def stress(lattice, lattice_inv, stress_3_3, rot):
    scaled_stress = np.dot(np.dot(lattice, stress_3_3), lattice.T)

    symmetrized_scaled_stress = np.zeros((3,3))
    for r in rot:
        symmetrized_scaled_stress += np.dot(np.dot(r.T,scaled_stress),r)
    symmetrized_scaled_stress /= len(rot)

    return np.dot(np.dot(lattice_inv, symmetrized_scaled_stress), lattice_inv.T)

# at = ase.io.read(sys.argv[1])
# 
# (rot, trans, symm_map) = prep(at)
# symm_for = force_symmetrize(at.get_cell(), at.get_reciprocal_cell().T, at.get_forces(), rot, trans, symm_map)
# symm_s = stress_symmetrize(at.get_cell(), at.get_reciprocal_cell().T, at.info['Virial'], rot)
# 
# print "forces before and after"
# for i_at in range(len(at)):
#     print at.get_forces()[i_at,:], symm_for[i_at,:]
# 
# print "stress before and after"
# print at.info['Virial']
# print symm_s
