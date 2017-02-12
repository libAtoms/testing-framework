import ase.io, os
from utilities import relax_config, model_test_root, run_root
import numpy as np

# the current 
import model 

def do_surface(test_dir):
    surf = ase.io.read(test_dir+"/surface.xyz", format="extxyz")

    # read bulk
    bulk_test_name=surf.info['bulk_test_name']
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

    # relax surface system
    tol = 1.0e-3
    surf = relax_config(surf, relax_pos=True, relax_cell=False, tol=tol, traj_file=None, config_label="surface", from_base_model=True, save_config=True)

    ase.io.write(os.path.join("..",run_root+"-relaxed.xyz"),  surf, format='extxyz')

    # check stoichiometry and number of bulk cell energies to subtract
    n_bulk_cells_0 = None
    surf_Zs = surf.get_atomic_numbers()
    bulk_Zs = bulk.get_atomic_numbers()
    for Z_i in set(bulk.get_atomic_numbers()):
        n_bulk_cells = sum(surf_Zs == Z_i)/sum(bulk_Zs == Z_i)
        if n_bulk_cells_0 is None:
            n_bulk_cells_0 = n_bulk_cells
        elif n_bulk_cells != n_bulk_cells_0:
            raise ValueError('surface has different stoichiometry from bulk')

    # calculate surface energy
    area = np.linalg.norm(np.cross(surf.get_cell()[0,:],surf.get_cell()[1,:]))
    return (bulk_model_test_root, (surf.get_potential_energy() - bulk.get_potential_energy()*n_bulk_cells)/(2.0*area) )
