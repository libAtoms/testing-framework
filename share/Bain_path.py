import sys
import ase.atoms
from utilities import get_relaxed_bulk, run_root, relax_config
import numpy as np
import spglib

def do_Bain_path(bulk_struct_test_fcc, bcc_fcc_steps=15, extra_deformation=0.3):

    c_a_step = (1.0 - 1.0/np.sqrt(2.0))/(bcc_fcc_steps-1)

    n_steps = bcc_fcc_steps
    # up
    c_a_t = 1.0
    n_extra_steps = int(np.ceil( (c_a_t*(1.0+extra_deformation)-c_a_t) / c_a_step ))
    n_steps += n_extra_steps
    c_a_max = 1.0 + c_a_step*n_extra_steps
    # down
    c_a_t = 1.0/np.sqrt(2)
    n_extra_steps = int(np.ceil( (c_a_t - c_a_t/(1.0+extra_deformation)) / c_a_step ))
    n_steps += n_extra_steps
    c_a_min = 1.0 / np.sqrt(2) - c_a_step*n_extra_steps

    properties = {}

    source_fcc = get_relaxed_bulk(bulk_struct_test_fcc)

    lattice, scaled_positions, numbers = spglib.standardize_cell(source_fcc)

    at0_fcc = ase.atoms.Atoms(cell=lattice, scaled_positions=scaled_positions, numbers=numbers, pbc=[True]*3)

    # assert normal cell vectors
    assert np.dot(at0_fcc.cell[0], at0_fcc.cell[1]) == 0
    assert np.dot(at0_fcc.cell[0], at0_fcc.cell[2]) == 0
    assert np.dot(at0_fcc.cell[1], at0_fcc.cell[2]) == 0
    # assert that cell[0] || \hat{x}, presumably a
    assert np.all(at0_fcc.cell[0][1:3] == 0)
    # assert that cell[0] || \hat{y}, presumably a
    assert np.all(at0_fcc.cell[1][[0,2]] == 0)
    # assert that cell[0] == cell[1]
    assert np.all(at0_fcc.cell[0][0] == at0_fcc.cell[1][1])
    # assert that cell[2] || \hat{z}, presumably c
    assert np.all(at0_fcc.cell[2][0:2] == 0)

    at = at0_fcc.copy()
    energies = []
    for c_a in np.linspace(c_a_max, c_a_min, n_steps):
        # fa*fa*fc = 1.0
        # fc*c / fa*a = c_a
        # fc = c_a * fa * a / c
        # 1.0 = fa * fa * (fa * c_a * a/c) = fa^3 c_a*a/c
        # fa = (c/(c_a*a))**(1/3)
        a = np.linalg.norm(at.cell[0])
        c = np.linalg.norm(at.cell[2])
        fa = (c/(c_a*a))**(1/3)
        fc = c_a * fa * a/c
        at.set_cell(np.diag((fa, fa, fc)) @ at.cell, True)
        at.info['c_a'] = c_a

        at = relax_config(at, relax_pos=False, relax_cell=True, hydrostatic_strain=True,
                          refine_symmetry_tol=1.0e-1, keep_symmetry=True,
                          config_label=f'Bain_path_{c_a:.3}', save_traj=True, from_base_model=True,
                          save_config=True)

        energies.append((c_a, at.get_potential_energy()/len(at), at.get_stress().tolist()))

    properties['Bain_path_E'] = energies

    return properties
