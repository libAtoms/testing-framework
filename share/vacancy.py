import ase.io, os.path
from utilities import relax_config, run_root
from ase.neighborlist import NeighborList
import numpy as np

def do_one_vacancy(relaxed_bulk, vac_i, relax_radial=0.0, relax_symm_break=0.0, nn_cutoff=0.0, tol=1.0e-3):
    vac = relaxed_bulk.copy()

    if relax_radial > 0:
        nl = NeighborList([nn_cutoff/2.0]*len(relaxed_bulk), self_interaction=False, bothways=True)
        nl.update(relaxed_bulk)
        indices, offsets = nl.get_neighbors(vac_i)
        offset_factor = relax_radial
        for i, offset in zip(indices, offsets):
           ri = vac.positions[vac_i] - (vac.positions[i] + np.dot(offset, vac.get_cell()))
           vac.positions[i] += offset_factor*ri
           offset_factor += relax_symm_break

    del vac[vac_i]
    label = "ind_%d_Z_%d" % (vac_i, relaxed_bulk.get_atomic_numbers()[vac_i])
    vac = relax_config(vac, relax_pos=True, relax_cell=False, tol=tol, traj_file=None, 
        config_label=label, from_base_model=True, save_config=True)

    ase.io.write(os.path.join("..",run_root+"-%s-relaxed.xyz" % label),  vac, format='extxyz')

    vac_pe = vac.get_potential_energy()
    if len(set(relaxed_bulk.get_atomic_numbers())) == 1:
        Ef0 = vac_pe - float(len(vac))/float(len(relaxed_bulk)) * relaxed_bulk.get_potential_energy()
    else:
        Ef0 =  vac_pe - relaxed_bulk.get_potential_energy()
    return ( label, run_root+"-%s-relaxed.xyz" % label, Ef0, relaxed_bulk.get_atomic_numbers()[vac_i] )

def do_all_vacancies(test_dir, nn_cutoff=0.0, tol=1.0e-3):
    print "doing do_all_vacancies"
    bulk = ase.io.read(os.path.join(test_dir,"bulk_supercell.xyz"), format="extxyz")
    print "got bulk ", len(bulk)

    tol = 1.0e-3
    relaxed_bulk = relax_config(bulk, relax_pos=True, relax_cell=True, tol=tol, 
        traj_file=None, config_label='bulk_supercell', from_base_model=True, save_config=True)
    relaxed_bulk_pe = relaxed_bulk.get_potential_energy()

    ase.io.write(os.path.join("..",run_root+"-relaxed-bulk.xyz"),  relaxed_bulk, format='extxyz')

    print "got relaxed_bulk ", relaxed_bulk.get_cell()

    properties = { "bulk_struct" : bulk.info["bulk_struct"], "bulk_E_per_atom" : relaxed_bulk_pe / len(relaxed_bulk), "defects" : {} }
    vacancy_list = [ int(i) for i in relaxed_bulk.info['vacancies'].split(',') ]

    for vac_i in vacancy_list:
        # maybe set up a system to read these from xyz file?
        relax_radial = 0.0
        relax_symm_break = 0.0
        (label, filename, Ef0, vac_Z) = do_one_vacancy(relaxed_bulk, vac_i, relax_radial, relax_symm_break, nn_cutoff, tol)
        if len(set(relaxed_bulk.get_atomic_numbers())) == 1:
            properties["defects"][label] = { 'Ef' : Ef0, 'filename' : filename, 'atom_ind' : vac_i, 'Z' : vac_Z }
        else:
            properties["defects"][label] = { 'Ef' : Ef0, 'filename' : filename, 'atom_ind' : vac_i, 'Z' : vac_Z, 'dmu' : [1, vac_Z] }

    return properties
