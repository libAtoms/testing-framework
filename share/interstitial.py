import ase.io, os.path
from ase import Atoms
from utilities import relax_config, run_root, rescale_to_relaxed_bulk, evaluate
from ase.neighborlist import NeighborList
from ase.constraints import FixedPlane
import numpy as np

def do_one_interstitial(bulk_supercell, bulk_supercell_pe, interstitial_Z, interstitial_pos, relax_radial=0.0, relax_symm_break=0.0, nn_cutoff=0.0, tol=1.0e-2):

    interst = bulk_supercell.copy()
    interst += Atoms(numbers=[interstitial_Z], positions=[interstitial_pos])
    interstitial_i = len(interst)-1

    if relax_radial != 0.0 or relax_symm_break != 0.0:
        nl = NeighborList([nn_cutoff/2.0]*len(bulk_supercell), self_interaction=False, bothways=True)
        nl.update(bulk_supercell)
        indices, offsets = nl.get_neighbors(interstitial_i)
        offset_factor = relax_radial
        for i, offset in zip(indices, offsets):
           ri = interst.positions[interstitial_i] - (interst.positions[i] + np.dot(offset, interst.get_cell()))
           interst.positions[i] += offset_factor*ri
           offset_factor += relax_symm_break

    if "interstitial_constraint" in bulk_supercell.info:
        (constr_type, constr_subtype) = bulk_supercell.info["interstitial_constraint"].split()[0:2]

        if constr_type == "plane":
            if constr_subtype == "atoms":
                indices = [int(i) for i in bulk_supercell.info["interstitial_constraint"].split()[2:]]
                if len(indices) != 3:
                    raise ValueError("number of indices not 3 for plane atoms '{}'".format(bulk_supercell.info["interstitial_constraint"]))
                p = interst.get_positions()
                constr_normal = np.cross(p[indices[0]]-p[indices[1]],p[indices[0]]-p[indices[2]])
            elif constr_subtype == "vector":
                constr_normal = np.array(bulk_supercell.info["interstitial_constraint"].split()[2:])
            else:
                raise ValueError("unknown interstitial constraint subtype for plane '{}'".format(bulk_supercell.info["interstitial_constraint"]))

            print "setting constraint FixedPlane with normal", constr_normal
            interst.set_constraint(FixedPlane(interstitial_i, constr_normal))
        else:
            raise ValueError("unknown interstitial constraint type '{}'".format(bulk_supercell.info["interstitial_constraint"]))

    evaluate(interst)
    unrelaxed_interstitial_pe = interst.get_potential_energy()

    label = "Z_%d" % (interstitial_Z)
    unrelaxed_filename=run_root+"-%s-unrelaxed.xyz" % label
    ase.io.write(os.path.join("..",unrelaxed_filename),  interst, format='extxyz')
    interst = relax_config(interst, relax_pos=True, relax_cell=False, tol=tol, traj_file=None, 
        config_label=label, from_base_model=True, save_config=True)

    relaxed_filename=run_root+"-%s-relaxed.xyz" % label
    ase.io.write(os.path.join("..",relaxed_filename),  interst, format='extxyz')

    interstitial_pe = interst.get_potential_energy()
    if len(set(bulk_supercell.get_atomic_numbers())) == 1:
        Ebulk = float(len(interst))/float(len(bulk_supercell)) * bulk_supercell_pe
    else:
        Ebulk = bulk_supercell_pe
    Ef0 = unrelaxed_interstitial_pe - Ebulk
    Ef = interstitial_pe - Ebulk
    print "got interstitial {} cell energy".format(label),interstitial_pe
    print "got bulk energy", Ebulk
    return ( label, unrelaxed_filename, Ef0, relaxed_filename, Ef, interstitial_i)


def do_interstitial(test_dir, nn_cutoff=0.0, tol=1.0e-2):
    print "doing do_interstitial"
    bulk_supercell = ase.io.read(os.path.join(test_dir,"bulk_supercell.xyz"), format="extxyz")
    print "got bulk_supercell ", len(bulk_supercell)

    bulk = rescale_to_relaxed_bulk(bulk_supercell)
    # relax bulk supercell positions in case it's only approximate (as it must be for different models), but stick 
    # to relaxed bulk's lattice constants as set by rescale_to_relaxed_bulk
    bulk_supercell = relax_config(bulk_supercell, relax_pos=True, relax_cell=False, tol=tol, traj_file=None, 
        config_label="relaxed_bulk", from_base_model=True, save_config=True)

    evaluate(bulk_supercell)
    bulk_supercell_pe = bulk_supercell.get_potential_energy()

    ase.io.write(os.path.join("..",run_root+"-rescaled-bulk.xyz"),  bulk_supercell, format='extxyz')

    print "got bulk primitive cell ", bulk.get_cell()
    print "got rescaled bulk_supercell cell ", bulk_supercell.get_cell()

    properties = { "bulk_struct_test" : bulk_supercell.info["bulk_struct_test"], "bulk_E_per_atom" : bulk_supercell_pe / len(bulk_supercell), "defects" : {} }

    try: # Cartesian 3-vector
        interstitial_pos = np.array([float(x) for x in bulk_supercell.info["interstitial_position"]])
        if len(interstitial_pos) != 3:
            raise ValueError("not a 3-vector")
    except:
        interstitial_pos_type = bulk_supercell.info["interstitial_position"].split()[0]
        if interstitial_pos_type == "mean":
            neighbor_indices = [int(i) for i in bulk_supercell.info["interstitial_position"].split()[1:]]
            if len(neighbor_indices) < 2:
                raise ValueError("interstitial position type mean, but {} < 2 indices".format(len(neighbor_indices)))
            interstitial_pos = np.sum(bulk_supercell.get_positions()[neighbor_indices],axis=0) / float(len(neighbor_indices))
        else:
            raise ValueError("Unknown interstitial position type in '"+bulk_supercell.info["interstitial_position"]+"'")

    if isinstance(bulk_supercell.info["Zs"],list):
        Z_list = bulk_supercell.info["Zs"]
    else:
        Z_list = [ bulk_supercell.info["Zs"] ]
    for interstitial_Z in Z_list:

        try:
            relax_radial = bulk_supercell.info['relax_radial_{}'.format(interstitial_Z)]
        except:
            relax_radial = 0.0
        try:
            relax_symm_break = bulk_supercell.info['relax_symm_break_{}'.format(interstitial_Z)]
        except:
            relax_symm_break = 0.0

        (label, unrelaxed_filename, Ef0, relaxed_filename, Ef, interstitial_i) = do_one_interstitial(bulk_supercell, bulk_supercell_pe, interstitial_Z, interstitial_pos, relax_radial, relax_symm_break, nn_cutoff, tol)

    properties["defects"][label] = { 'Ef0' : Ef0, 'Ef' : Ef, 'unrelaxed_filename' : unrelaxed_filename, 'relaxed_filename' : relaxed_filename, 'atom_ind' : interstitial_i, 'Z' : interstitial_Z }
    if len(set(bulk_supercell.get_atomic_numbers())) != 1:
        properties["defects"][label]['dmu'] = [-1, interstitial_Z]

    return properties
