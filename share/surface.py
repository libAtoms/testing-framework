import ase.io, os
from utilities import relax_config, model_test_root, run_root
import numpy as np

# the current 
import model 

def do_symmetric_surface(test_dir):
    surf = ase.io.read(test_dir+"/surface.xyz", format="extxyz")

    # read bulk
    bulk_test_name=surf.info['bulk_struct']
    bulk_model_test_root = model_test_root(u_test_name=bulk_test_name)

    bulk = ase.io.read('../%s-relaxed.xyz' % bulk_model_test_root, format='extxyz')

    # rescale surface cell
    surf_a1_lattice = surf.info['surf_a1_in_lattice_coords']
    bulk_lattice = bulk.get_cell()
    surf_a1_in_bulk = np.dot(surf_a1_lattice,bulk_lattice)
    cell_ratio = np.linalg.norm(surf_a1_in_bulk) / np.linalg.norm(surf.get_cell()[0,:])

    if 'surf_a2_in_lattice_coords' in surf.info:
        raise ValueError('anisotropic rescaling of surface cell not implemented')
    if 'surf_a3_in_lattice_coords' in surf.info:
        raise ValueError('anisotropic rescaling of surface cell not implemented')

    surf.set_cell(surf.get_cell()*cell_ratio, scale_atoms=True)

    print "got relaxed bulk cell ", bulk.get_cell()
    print "got rescaled surf cell ", surf.get_cell()

    # relax surface system
    tol = 1.0e-3
    surf = relax_config(surf, relax_pos=True, relax_cell=False, tol=tol, traj_file=None, config_label="surface", from_base_model=True, save_config=True)

    ase.io.write(os.path.join("..",run_root+"-relaxed.xyz"),  surf, format='extxyz')

    # check stoichiometry and number of bulk cell energies to subtract
    surf_Zs = surf.get_atomic_numbers()
    bulk_Zs = bulk.get_atomic_numbers()
    Z0 = bulk_Zs[0]
    n_bulk_cells = float(sum(surf_Zs == Z0))/float(sum(bulk_Zs == Z0))
    if len(set(bulk_Zs)) == 1:
        n_dmu = None
    else:
        n_dmu = {}
        for Z in set(bulk_Zs):
            n_dmu[Z] = n_bulk_cells*sum(bulk_Zs == Z) - sum(surf_Zs == Z)

    # calculate surface energy
    area = np.linalg.norm(np.cross(surf.get_cell()[0,:],surf.get_cell()[1,:]))
    return { "bulk_struct" : bulk_test_name,  "Ef" : (surf.get_potential_energy() - bulk.get_potential_energy()*n_bulk_cells)/(2.0*area), "dmu" : n_dmu }
