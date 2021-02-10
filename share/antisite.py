import ase.io, os.path
from utilities import relax_config, run_root, rescale_to_relaxed_bulk, evaluate
import ase.build
import numpy as np
import spglib

def do_one_antisite_pair(bulk_supercell, bulk_supercell_pe, i1, i2, tol=1.0e-2):
    assert bulk_supercell.numbers[i1] != bulk_supercell.numbers[i2]

    # do unrelaxed (without perturbations)
    antisite = bulk_supercell.copy()
    Z1, Z2 = antisite.numbers[[i1, i2]]
    antisite.numbers[i1] = Z2
    antisite.numbers[i2] = Z1

    label = "ind_%d_Z_%d_ind_%d_Z_%d" % (i1, Z1, i2, Z2)
    unrelaxed_filename=run_root+"-%s-unrelaxed.xyz" % label
    ase.io.write(os.path.join("..", unrelaxed_filename), antisite, format='extxyz')
    evaluate(antisite)
    unrelaxed_antisite_pe = antisite.get_potential_energy()

    antisite = relax_config(antisite, relax_pos=True, relax_cell=False, tol=tol, save_traj=True,
        config_label=label, from_base_model=True, save_config=True, try_restart=True)
    relaxed_filename=run_root+"-%s-relaxed.xyz" % label
    ase.io.write(os.path.join("..", relaxed_filename), antisite, format='extxyz')

    # already has calculator from relax_configs
    antisite_pe = antisite.get_potential_energy()
    Ef0 = unrelaxed_antisite_pe - bulk_supercell_pe
    Ef = antisite_pe - bulk_supercell_pe
    print("got antisite", label, "cell energy", antisite_pe, "n_atoms", len(antisite))
    print("got bulk energy", bulk_supercell_pe)
    return ( label, unrelaxed_filename, Ef0, relaxed_filename, Ef, Z1, Z2 )

def do_farthest_inequiv_pairs(test_dir, tol=1.0e-2):
    print("doing do_farthest_antisite_pairs")
    bulk_supercell = ase.io.read(os.path.join(test_dir, "bulk_supercell.xyz"), format="extxyz")
    print("got bulk_supercell ", len(bulk_supercell))

    bulk = rescale_to_relaxed_bulk(bulk_supercell)
    # relax bulk supercell positions in case it's only approximate (as it must be for different models), but stick
    # to relaxed bulk's lattice constants as set by rescale_to_relaxed_bulk
    bulk_supercell = relax_config(bulk_supercell, relax_pos=True, relax_cell=False, tol=tol, save_traj=True,
        config_label="rescaled_bulk", from_base_model=True, save_config=True)

    ase.io.write(os.path.join("..", run_root+"-rescaled-bulk.xyz"),  bulk_supercell, format='extxyz')

    print("got bulk primitive cell ", bulk.get_cell())
    print("got rescaled bulk_supercell cell ", bulk_supercell.get_cell())

    if 'arb_supercell' in bulk_supercell.info:
        print("making bulk supercell from", bulk_supercell.info['arb_supercell'].reshape((3,3)))
        bulk_supersupercell = ase.build.make_supercell(bulk_supercell, bulk_supercell.info['arb_supercell'].reshape((3,3)))
        print("got supersupercell with ", len(bulk_supersupercell), "atoms, cell\n", bulk_supersupercell.get_cell())

        bulk_supersupercell.info.update(bulk_supercell.info)
        bulk_supercell = bulk_supersupercell

    sym_data = spglib.get_symmetry_dataset(bulk_supercell, symprec=0.01)
    equiv_at = set([tuple(iZ) for iZ in zip(sym_data["equivalent_atoms"], bulk_supercell.numbers)])

    # print("equiv_at", equiv_at)

    antisite_list = []
    for i1, Z1 in equiv_at:
        for i2_proto, Z2 in equiv_at:
            if Z1 <= Z2:
                continue
            # print("check i", i1, i2_proto, "Z", Z1, Z2)

            i2s = np.where(sym_data["equivalent_atoms"] == i2_proto)[0]
            i2_dists = bulk_supercell.get_distances(i1, i2s, mic=True)
            farthest_ind = np.argmax(i2_dists)
            i2 = i2s[farthest_ind]

            antisite_list.append((i1, i2))

    # print("antisite_list", antisite_list) ##

    evaluate(bulk_supercell)
    bulk_supercell_pe = bulk_supercell.get_potential_energy()
    properties = { "bulk_struct_test" : bulk_supercell.info["bulk_struct_test"], "bulk_E_per_atom" : bulk_supercell_pe / len(bulk_supercell), "defects" : {} }

    for i1, i2 in antisite_list:
        (label, unrelaxed_filename, Ef0, relaxed_filename, Ef, Z1, Z2) = do_one_antisite_pair(bulk_supercell, bulk_supercell_pe, i1, i2, tol)

        properties["defects"][label] = { 'Ef0' : Ef0, 'Ef' : Ef, 'unrelaxed_filename' : unrelaxed_filename, 'relaxed_filename' : relaxed_filename,
            'atom_inds' : (int(i1), int(i2)), 'Zs' : (int(Z1), int(Z2)) }

    print("returning properties", properties)
    return properties
