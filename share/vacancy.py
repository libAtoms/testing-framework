import ase.io, os.path
from utilities import relax_config, run_root

def do_one_vacancy(relaxed_bulk, vac_i, tol=1.0e-3):
    vac = relaxed_bulk.copy()
    del vac[vac_i]
    label = "ind_%d_Z_%d" % (vac_i, relaxed_bulk.get_atomic_numbers()[vac_i])
    vac = relax_config(vac, relax_pos=True, relax_cell=False, tol=tol, traj_file=None, 
        label=label, from_base_model=True, save_config=True)

    ase.io.write(os.path.join("..",run_root+"-relaxed.xyz"),  vac, format='extxyz')

    vac_pe = vac.get_potential_energy()
    if len(set(vac.get_atomic_numbers())) == 1:
        return [label, vac_pe - float(len(vac))/float(len(relaxed_bulk)) * relaxed_bulk.get_potential_energy() ]
    else:
        return [label, vac_pe - relaxed_bulk.get_potential_energy(), "+mu_%d" % relaxed_bulk.get_atomic_numbers()[vac_i] ]

def do_all_vacancies(test_dir, tol=1.0e-3):
    bulk = ase.io.read(os.path.join(test_dir,"bulk_supercell.xyz"), format="extxyz")

    tol = 1.0e-3
    bulk = relax_config(bulk, relax_pos=True, relax_cell=True, tol=tol, 
        traj_file=None, label='bulk_supercell', from_base_model=True, save_config=True)

    properties={}
    try:
        vacancy_list = [ int(i) for i in bulk.info['vacancies'].split(',') ]
    except:
        vacancy_list = [ bulk.info['vacancies'] ]
    for vac_i in vacancy_list:
        vac_E = do_one_vacancy(bulk, vac_i, tol)
        if len(vac_E) == 2:
            properties[vac_E[0]] = vac_E[1]
        elif len(vac_E) == 3:
            properties[vac_E[0]] = [vac_E[1],vac_E[2]]
        else:
            raise ValueError('unknown number of things returned from do_vacancy')

    return properties
