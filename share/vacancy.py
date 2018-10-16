import ase.io, os.path
from utilities import relax_config, run_root, rescale_to_relaxed_bulk, evaluate
from ase.neighborlist import NeighborList
import numpy as np

def do_one_vacancy(bulk_supercell, bulk_supercell_pe, vac_i, relax_radial=0.0, relax_symm_break=0.0, nn_cutoff=0.0, tol=1.0e-3):

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
    vac = relax_config(vac, relax_pos=True, relax_cell=False, tol=tol, traj_file=None, 
        config_label=label, from_base_model=True, save_config=True)

    ase.io.write(os.path.join("..",run_root+"-%s-relaxed.xyz" % label),  vac, format='extxyz')

    vac_pe = vac.get_potential_energy()
    if len(set(bulk_supercell.get_atomic_numbers())) == 1:
        Ef0 = vac_pe - float(len(vac))/float(len(bulk_supercell)) * bulk_supercell_pe
    else:
        Ef0 =  vac_pe - bulk_supercell_pe
    return ( label, run_root+"-%s-relaxed.xyz" % label, Ef0, bulk_supercell.get_atomic_numbers()[vac_i] )

def do_all_vacancies(test_dir, nn_cutoff=0.0, tol=1.0e-3):
    print "doing do_all_vacancies"
    bulk_supercell = ase.io.read(os.path.join(test_dir,"bulk_supercell.xyz"), format="extxyz")
    print "got bulk_supercell ", len(bulk_supercell)

    bulk = rescale_to_relaxed_bulk(bulk_supercell)

    tol = 1.0e-3
    evaluate(bulk_supercell)
    bulk_supercell_pe = bulk_supercell.get_potential_energy()

    ase.io.write(os.path.join("..",run_root+"-rescaled-bulk.xyz"),  bulk_supercell, format='extxyz')

    print "got bulk primitive cell ", bulk.get_cell()
    print "got rescaled bulk_supercell cell ", bulk_supercell.get_cell()

    properties = { "bulk_struct_test" : bulk_supercell.info["bulk_struct_test"], "bulk_E_per_atom" : bulk_supercell_pe / len(bulk_supercell), "defects" : {} }
    try:
        vacancy_list = [ int(i) for i in bulk_supercell.info['vacancies'] ]
    except:
        vacancy_list = [ int(bulk_supercell.info['vacancies']) ]

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
        (label, filename, Ef0, vac_Z) = do_one_vacancy(bulk_supercell, bulk_supercell_pe, vac_i, relax_radial, relax_symm_break, nn_cutoff, tol)
        if len(set(bulk_supercell.get_atomic_numbers())) == 1:
            properties["defects"][label] = { 'Ef' : Ef0, 'filename' : filename, 'atom_ind' : vac_i, 'Z' : vac_Z }
        else:
            properties["defects"][label] = { 'Ef' : Ef0, 'filename' : filename, 'atom_ind' : vac_i, 'Z' : vac_Z, 'dmu' : [1, vac_Z] }

    return properties
