from utilities import *
import numpy as np
import os
import phonopy
from phonopy.file_IO import parse_FORCE_CONSTANTS, write_FORCE_CONSTANTS
import ase.units

def do_phonons(bulk_struct_tests, n_supercell, n_dos_mesh=0, dx=0.01):
    properties = {}
    for bulk_struct_test in bulk_struct_tests:
        at0 = get_relaxed_bulk(bulk_struct_test)

        # magnetic moments could change the symmetry, ignored here for now
        phonopy_atoms = phonopy.structure.atoms.PhonopyAtoms( symbols=at0.get_chemical_symbols(), 
                                      scaled_positions=at0.get_scaled_positions(),
                                      masses=at0.get_masses(), cell=at0.get_cell() )

        phonons = phonopy.Phonopy( phonopy_atoms, np.diag([n_supercell]*3), factor=phonopy.units.VaspToTHz )
        phonons.generate_displacements(distance=dx)

        # convert from chosen Phonopy units (THz) to cm^-1
        THz_per_invcm = ase.units._c * 1.0e-10

        ####################################################################################################
        # if args.SETUP:

        displ_supercells = phonons.get_supercells_with_displacements()

        at0 *= n_supercell
        # reorder at0 to match order from phonopy
        at = ase.Atoms(pbc=True, cell=displ_supercells[0].get_cell(), positions=displ_supercells[0].get_positions(), numbers=displ_supercells[0].get_atomic_numbers())
        matched_pos = np.zeros(at.positions.shape)
        mapping = [-1]*len(at)
        at0_scaled_pos = at0.get_scaled_positions()
        at_scaled_pos = at.get_scaled_positions()
        for at0_i in range(len(at)):
            scaled_dists = at0_scaled_pos[at0_i] - at_scaled_pos
            scaled_dists -= np.round(scaled_dists)
            closest_i = np.argmin(np.linalg.norm(scaled_dists, axis=1))
            matched_pos[at0_i] = at.positions[closest_i]
            mapping[closest_i] = at0_i
        if -1 in mapping:
            raise RuntimeError("Failed to map orig and displaced atom positions")
        at0 = at0[mapping]

        # ase.io.write("UNDISPL.{}".format(file_label), at0)

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

        sys.stderr.write("Calculating for {} displacements\n".format(len(at_sets)))

        evaluate(at0)
        f0 = at0.get_forces()
        # np.savetxt("FORCES.UNDISPL.{}".format(file_label), f0)
        p0 = at0.get_positions()
        c0 = at0.get_cell()

        for (displ_i, at) in enumerate(at_sets):
            at0.set_positions(at.positions)
            at0.set_cell(at.get_cell())
            evaluate(at0)
            f = at0.get_forces()
            all_forces.append(f)
            # np.savetxt("FORCES.DISPL_{}.{}".format(displ_i, file_label), all_forces[-1])

        ####################################################################################################
        # if args.PROCESS:
        Nat = all_forces[0].shape[0]
        for f in all_forces:
            f -= f0
            f -= np.outer(np.ones(Nat), np.sum(f, axis=0)/Nat)

        phonons.produce_force_constants(forces=all_forces, calculate_full_force_constants=True )
        phonons.symmetrize_force_constants()

        # CsG says a bug in phonopy requires writing and then re-reading the force constants
        write_FORCE_CONSTANTS(phonons.get_force_constants(), filename="FORCE_CONSTANTS")

        ####################################################################################################
        # if args.ANALYZE:
        # CsG says an apparent bug in phonopy requires that phonons be written and then read back in, so don't just use the ones calculated in prev stage directly
        fc = parse_FORCE_CONSTANTS(filename="FORCE_CONSTANTS")
        phonons.set_force_constants(fc)

        # do something with force constants
        properties[bulk_struct_test] = {}
        if n_dos_mesh > 0:
            phonons.set_mesh( [n_dos_mesh]*3 )

            phonons.set_total_DOS()
            frequencies = phonons.get_total_dos_dict()['frequency_points']
            PHdos = phonons.get_total_dos_dict()['total_dos']

            #contains by columns the frequencies in cm^{-1} and the vDOS
            #the vDOS ins in units of "number of states/(unit cell x frequency[cm^{-1}])" 
            # i.e. if you integrate the vDOS throughout frequency, it will be 3N, where N is the number of atoms in the unit cell
            # np.savetxt( "DOS.{}{}.data".format(file_label,calculator_str), np.transpose([frequencies/THz_per_invcm,PHdos*THz_per_invcm]))
            properties[bulk_struct_test]['DOS'] = [list(frequencies/THz_per_invcm), list(PHdos*THz_per_invcm)]

        ####################################################################################################

        return properties
