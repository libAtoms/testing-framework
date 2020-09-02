import sys
from utilities import *
import numpy as np
import phonopy
import ase.units

def do_phonons(bulk_struct_tests, n_supercell, band_paths=None, dx=0.01):
    if band_paths is not None and len(band_paths) != len(bulk_struct_tests):
        raise RuntimeError("got {} bulk structs but different {} band paths".format(len(bulk_struct_tests), len(band_paths)))
    properties = {}
    if type(n_supercell) != list:
        n_supercell_list = [ n_supercell for i in bulk_struct_tests ]
    else:
        n_supercell_list = n_supercell
    
    for bulk_i, bulk_struct_test in enumerate(bulk_struct_tests):
        at0 = get_relaxed_bulk(bulk_struct_test)

        # magnetic moments could change the symmetry, ignored here for now
        phonopy_atoms = phonopy.structure.atoms.PhonopyAtoms( symbols=at0.get_chemical_symbols(), 
                                      scaled_positions=at0.get_scaled_positions(),
                                      masses=at0.get_masses(), cell=at0.get_cell() )

        if np.size(n_supercell_list[bulk_i]) == 1:
            phonons = phonopy.Phonopy( phonopy_atoms, np.diag([n_supercell_list[bulk_i]]*3), factor=phonopy.units.VaspToTHz )
        else:
            phonons = phonopy.Phonopy( phonopy_atoms, np.diag(n_supercell_list[bulk_i]), factor=phonopy.units.VaspToTHz )	
        phonons.generate_displacements(distance=dx)

        # convert from chosen Phonopy units (THz) to cm^-1
        THz_per_invcm = ase.units._c * 1.0e-10

        ####################################################################################################
        # if args.SETUP:
        sys.stderr.write("SETUP\n")

        displ_supercells = phonons.get_supercells_with_displacements()

        at0_sc = at0 * n_supercell_list[bulk_i]
        # reorder at0_sc to match order from phonopy
        at = ase.Atoms(pbc=True, cell=displ_supercells[0].get_cell(), positions=displ_supercells[0].get_positions(), numbers=displ_supercells[0].get_atomic_numbers())
        matched_pos = np.zeros(at.positions.shape)
        mapping = [-1]*len(at)
        at0_sc_scaled_pos = at0_sc.get_scaled_positions()
        at_scaled_pos = at.get_scaled_positions()
        for at0_sc_i in range(len(at)):
            scaled_dists = at0_sc_scaled_pos[at0_sc_i] - at_scaled_pos
            scaled_dists -= np.round(scaled_dists)
            closest_i = np.argmin(np.linalg.norm(scaled_dists, axis=1))
            matched_pos[at0_sc_i] = at.positions[closest_i]
            mapping[closest_i] = at0_sc_i
        if -1 in mapping:
            raise RuntimeError("Failed to map orig and displaced atom positions")
        at0_sc = at0_sc[mapping]

        # ase.io.write("UNDISPL.{}".format(file_label), at0_sc)

        # create displaced cells
        sys.stderr.write("Creating {} displacements\n".format(len(displ_supercells)))
        at_sets = []
        for (displ_i, displ_config) in enumerate(displ_supercells):
            at = ase.Atoms(pbc=True, cell=displ_config.get_cell(), positions=displ_config.get_positions(), numbers=displ_config.get_atomic_numbers())
            at_sets.append(at)
            # ase.io.write("DISPL_{}.{}".format(displ_i, file_label), at)

        ####################################################################################################
        all_forces = []
        # if args.CALCULATE:
        sys.stderr.write("CALCULATE\n")

        sys.stderr.write("Calculating for {} displacements\n".format(len(at_sets)))

        evaluate(at0_sc)
        f0 = at0_sc.get_forces()
        # np.savetxt("FORCES.UNDISPL.{}".format(file_label), f0)

        for (displ_i, at) in enumerate(at_sets):
            at0_sc.set_positions(at.positions)
            at0_sc.set_cell(at.get_cell())
            sys.stderr.write("start evaluate displacement {}\n". format(displ_i))
            evaluate(at0_sc)
            f = at0_sc.get_forces()
            all_forces.append(f)
            # np.savetxt("FORCES.DISPL_{}.{}".format(displ_i, file_label), all_forces[-1])

        ####################################################################################################
        # if args.PROCESS:
        sys.stderr.write("PROCESS\n")

        Nat = all_forces[0].shape[0]
        for f in all_forces:
            f -= f0
            f -= np.outer(np.ones(Nat), np.sum(f, axis=0)/Nat)
        all_forces = np.asarray(all_forces)

        properties[bulk_struct_test] = { 'dx' : dx, 'all_forces' : all_forces.tolist(), 'symb' : at0.get_chemical_symbols(), 'scaled_pos' : at0.get_scaled_positions().tolist(), 'm' : at0.get_masses().tolist(), 'c' : at0.get_cell().tolist(), 'n_cell' : np.diag([n_supercell_list[bulk_i]]*3).tolist(), 'unit_factor' : phonopy.units.VaspToTHz, 'band_path' : band_paths[bulk_i] }

        ####################################################################################################

    return properties
