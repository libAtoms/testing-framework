import ase.io, os.path
from utilities import relax_config, run_root, rescale_to_relaxed_bulk, evaluate
from ase.neighborlist import NeighborList
import ase.build
import numpy as np
import spglib

def do_one_vacancy(bulk_supercell, bulk_supercell_pe, vac_i, relax_radial=0.0, relax_symm_break=0.0, nn_cutoff=0.0, tol=1.0e-2):

    vac = bulk_supercell.copy()

    if relax_radial != 0.0 or relax_symm_break != 0.0:
        nl = NeighborList([nn_cutoff/2.0]*len(bulk_supercell), self_interaction=False, bothways=True)
        nl.update(bulk_supercell)
        indices, offsets = nl.get_neighbors(vac_i)
        offset_factor = relax_radial
        for i, offset in zip(indices, offsets):
           ri = vac.positions[vac_i] - (vac.positions[i] + np.dot(offset, vac.get_cell()))
           vac.positions[i] += offset_factor*ri
           offset_factor += relax_symm_break

    del vac[vac_i]

    label = "ind_%d_Z_%d" % (vac_i, bulk_supercell.get_atomic_numbers()[vac_i])
    unrelaxed_filename=run_root+"-%s-unrelaxed.xyz" % label
    ase.io.write(os.path.join("..",unrelaxed_filename), vac, format='extxyz')
    evaluate(vac)
    unrelaxed_vac_pe = vac.get_potential_energy()

    vac = relax_config(vac, relax_pos=True, relax_cell=False, tol=tol, traj_file=None, 
        config_label=label, from_base_model=True, save_config=True)
    relaxed_filename=run_root+"-%s-relaxed.xyz" % label
    ase.io.write(os.path.join("..",relaxed_filename), vac, format='extxyz')

    # already has calculator from relax_configs
    vac_pe = vac.get_potential_energy()
    if len(set(bulk_supercell.get_atomic_numbers())) == 1:
        Ebulk = float(len(vac))/float(len(bulk_supercell)) * bulk_supercell_pe
    else:
        Ebulk = bulk_supercell_pe
    Ef0 = unrelaxed_vac_pe - Ebulk
    Ef = vac_pe - Ebulk
    print "got vacancy",label,"cell energy",vac_pe,"n_atoms",len(vac)
    print "got bulk energy", Ebulk," (scaled to (N-1)/N if single component)"
    return ( label, unrelaxed_filename, Ef0, relaxed_filename, Ef, bulk_supercell.get_atomic_numbers()[vac_i] )

def do_all_vacancies(test_dir, nn_cutoff=0.0, tol=1.0e-2):
    print "doing do_all_vacancies"
    bulk_supercell = ase.io.read(os.path.join(test_dir,"bulk_supercell.xyz"), format="extxyz")
    print "got bulk_supercell ", len(bulk_supercell)

    bulk = rescale_to_relaxed_bulk(bulk_supercell)
    # relax bulk supercell positions in case it's only approximate (as it must be for different models), but stick 
    # to relaxed bulk's lattice constants as set by rescale_to_relaxed_bulk
    bulk_supercell = relax_config(bulk_supercell, relax_pos=True, relax_cell=False, tol=tol, traj_file=None, 
        config_label="rescaled_bulk", from_base_model=True, save_config=True)

    ase.io.write(os.path.join("..",run_root+"-rescaled-bulk.xyz"),  bulk_supercell, format='extxyz')

    print "got bulk primitive cell ", bulk.get_cell()
    print "got rescaled bulk_supercell cell ", bulk_supercell.get_cell()

    if bulk_supercell.info['vacancies'] == "inequivalent":
        sym_data = spglib.get_symmetry_dataset(bulk_supercell, symprec=0.01)
        prim_vacancy_list = np.unique(sym_data["equivalent_atoms"])
        print "orig cell vacancy_list", prim_vacancy_list
        if 'arb_supercell' in bulk_supercell.info:
            bulk_supersupercell = ase.build.cut(bulk_supercell, bulk_supercell.info['arb_supercell'].T[0], 
                bulk_supercell.info['arb_supercell'].T[1], bulk_supercell.info['arb_supercell'].T[2])
            print "get supersupercell with ",len(bulk_supersupercell),"atoms, cell\n",bulk_supersupercell.get_cell()
            vacancy_list = []
            for i in prim_vacancy_list:
                p = bulk_supercell.get_positions()[i]
                dv = bulk_supersupercell.get_positions() - p
                dv_scaled = np.dot(dv, bulk_supersupercell.get_reciprocal_cell().T)
                dv -= np.dot(np.round(dv_scaled), bulk_supersupercell.get_cell())
                i_closest = np.argmin(np.linalg.norm(dv, axis=1))
                print "found closest in new cell", i_closest, "distance in orig cell lattice coords", np.dot((bulk_supersupercell.get_positions()[i_closest]-p), 
                    bulk_supercell.get_reciprocal_cell().T)
                vacancy_list.append(i_closest)
            bulk_supercell = bulk_supersupercell
        else:
            vacancy_list = prim_vacancy_list
        print "final vacancy_list", vacancy_list
    else:
        try:
            vacancy_list = [ int(i) for i in bulk_supercell.info['vacancies'] ]
        except:
            vacancy_list = [ int(bulk_supercell.info['vacancies']) ]

    evaluate(bulk_supercell)
    bulk_supercell_pe = bulk_supercell.get_potential_energy()
    properties = { "bulk_struct_test" : bulk_supercell.info["bulk_struct_test"], "bulk_E_per_atom" : bulk_supercell_pe / len(bulk_supercell), "defects" : {} }

    for vac_i in vacancy_list:
        # maybe set up a system to read these from xyz file?
        try:
            relax_radial = bulk_supercell.info['relax_radial_{}'.format(vac_i)]
        except:
            relax_radial = 0.0
        try:
            relax_symm_break = bulk_supercell.info['relax_symm_break_{}'.format(vac_i)]
        except:
            relax_symm_break = 0.0
        (label, unrelaxed_filename, Ef0, relaxed_filename, Ef, vac_Z) = do_one_vacancy(bulk_supercell, bulk_supercell_pe, vac_i, relax_radial, relax_symm_break, nn_cutoff, tol)

        properties["defects"][label] = { 'Ef0' : Ef0, 'Ef' : Ef, 'unrelaxed_filename' : unrelaxed_filename,'relaxed_filename' : relaxed_filename, 'atom_ind' : vac_i, 'Z' : vac_Z }
        if len(set(bulk_supercell.get_atomic_numbers())) > 1:
            properties["defects"][label]["dmu"] = [1, vac_Z]

    return properties
