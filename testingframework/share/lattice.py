import numpy as np
from testingframework.share.utilities  import relax_config, run_root
import ase.io, sys, os.path
from ase.optimize.precon import PreconLBFGS
import matscipy.elasticity
from ase.units import GPa
# Calculate B for Hexagonal, Tetragonal & Trigonal
# DOI: 10.1103/PhysRevB.77.104118 Eq. (29)
def HTT_B(c11, c33, c12, c13):
    numerator   = c33 * (c11 + c12) - 2*c13**2
    denominator = c11 + c12 + 2*c33 - 4*c13

    return numerator / denominator

# Mater Trans v 53 p 1247 (2012), Eq. 4-10
def VRH_B(c11, c33, c12, c13, c44, c66):
    Bv = (1.0/9.0)*(2.0*(c11+c12) + c33 + 4.0*c13)

    M = c11 + c12 + 2.0*c33 - 4.0*c13
    Csq = (c11+c12)*c33 - 2.0*c13**2
    Br = Csq/M

    return (Bv+Br)/2.0

def calc_E_vs_V(bulk, dV=0.025, n_steps=(-10,10), tol=1e-2, method='lbfgs'): # hack tol to deal with Te C2/m
   import model

   V0 = bulk.get_volume()
   dV *= V0
   E_vs_V=[]

   scaled_bulk = bulk.copy()
   for i in range(0, n_steps[0]-1, -1):
      V_cur = scaled_bulk.get_volume()
      scaled_bulk.set_cell(scaled_bulk.get_cell()*((V0+i*dV)/V_cur)**(1.0/3.0), scale_atoms=True)
      ase.io.write(sys.stdout, scaled_bulk, format='extxyz')
      print("trying to relax i",i)
      try:
          if hasattr(model, "fix_cell_dependence"):
               model.fix_cell_dependence(scaled_bulk)
          ase.io.write(run_root+"-E_vs_V_%03d-unrelaxed.xyz" % i,  scaled_bulk, format='extxyz')
          scaled_bulk = relax_config(scaled_bulk, relax_pos=True, relax_cell=True, tol=tol, max_steps=200, save_traj=True, constant_volume=True, method=method,
              refine_symmetry_tol=1.0e-1, keep_symmetry=True, config_label="E_vs_V_%03d" % i, from_base_model=True, save_config=True)
      except Exception as e:
          print("WARNING: failed config in calc_E_vs_V", str(e))
          sys.exit(1) #### NB
          break
      ase.io.write(sys.stdout, scaled_bulk, format='extxyz')
      E_vs_V.insert(0, (scaled_bulk.get_volume()/len(scaled_bulk), scaled_bulk.get_potential_energy()/len(bulk), list(scaled_bulk.get_stress())) )

   scaled_bulk = bulk.copy()
   for i in range(1,n_steps[1]+1):
      V_cur = scaled_bulk.get_volume()
      scaled_bulk.set_cell(scaled_bulk.get_cell()*((V0+i*dV)/V_cur)**(1.0/3.0), scale_atoms=True)
      ase.io.write(sys.stdout, scaled_bulk, format='extxyz')
      print("trying to relax i",i)
      try:
          if hasattr(model, "fix_cell_dependence"):
               model.fix_cell_dependence(scaled_bulk)
          ase.io.write(run_root+"-E_vs_V_%02d-unrelaxed.xyz" % i,  scaled_bulk, format='extxyz')
          scaled_bulk = relax_config(scaled_bulk, relax_pos=True, relax_cell=True, tol=tol, max_steps=200, save_traj=True, constant_volume=True, method=method,
              refine_symmetry_tol=1.0e-1, keep_symmetry=True, config_label="E_vs_V_%02d" % i, from_base_model=True, save_config=True)
      except Exception as e:
          print("failed", str(e))
          break
      ase.io.write(sys.stdout, scaled_bulk, format='extxyz')
      E_vs_V.append( (scaled_bulk.get_volume()/len(scaled_bulk), scaled_bulk.get_potential_energy()/len(bulk), list(scaled_bulk.get_stress())) )

   if hasattr(model, "fix_cell_dependence"):
       model.fix_cell_dependence()

   return E_vs_V


def do_lattice(test_dir, lattice_type, dV=0.025, n_steps=(-10,10), tol=1.0e-2, method='lbfgs', applied_P=0.0):

   import model
   bulk = ase.io.read(test_dir+"/bulk.xyz", format="extxyz")

   results_dict = {}

   print("relax bulk")
   # relax the initial unit cell and atomic positions
   (orig_cell, new_cell) = (None, None)
   while new_cell is None or np.max(np.abs(np.dot(np.linalg.inv(new_cell),orig_cell) - np.eye(3))) > 0.05:
       if hasattr(model, "fix_cell_dependence"):
           model.fix_cell_dependence(bulk)
       orig_cell = bulk.get_cell()
       bulk = relax_config(bulk, relax_pos=True, relax_cell=True, tol=tol, save_traj=True, method=method,
                           refine_symmetry_tol=1.0e-2, keep_symmetry=True, config_label="bulk", from_base_model=True, save_config=True, applied_P=applied_P)
       new_cell = bulk.get_cell()
       if hasattr(model, "fix_cell_dependence"):
           model.fix_cell_dependence()
       else:
           break

   print("final relaxed bulk")
   ase.io.write(sys.stdout, bulk, format='extxyz')
   ase.io.write(os.path.join("..",run_root+"-relaxed.xyz"),  bulk, format='extxyz')

   print("calculating E vs. V")
   E_vs_V = calc_E_vs_V(bulk, dV=dV, n_steps=n_steps, tol=tol)
   results_dict.update({ 'E_vs_V' : E_vs_V })

   print("calculating elastic constants")

   if hasattr(model, "fix_cell_dependence"):
       model.fix_cell_dependence(bulk)

   opt = lambda atoms, **kwargs: PreconLBFGS(atoms, **kwargs)
   if lattice_type == 'cubic':
       elastic_consts = matscipy.elasticity.fit_elastic_constants(bulk, symmetry='cubic', optimizer=opt, logfile=sys.stdout)
       c11 = elastic_consts[0][0,0]/GPa
       c12 = elastic_consts[0][0,1]/GPa
       c44 = elastic_consts[0][3,3]/GPa
       results_dict.update({'c11' : c11, 'c12': c12, 'c44' : c44, 'B' : (c11+2.0*c12)/3.0})
   elif lattice_type == 'orthorhombic':
       elastic_consts = matscipy.elasticity.fit_elastic_constants(bulk, optimizer=opt, logfile=sys.stdout)
       c11 = elastic_consts[0][0,0]/GPa
       c22 = elastic_consts[0][1,1]/GPa
       c33 = elastic_consts[0][2,2]/GPa
       c12 = elastic_consts[0][0,1]/GPa
       c13 = elastic_consts[0][0,2]/GPa
       c23 = elastic_consts[0][1,2]/GPa
       c44 = elastic_consts[0][3,3]/GPa
       c55 = elastic_consts[0][4,4]/GPa
       c66 = elastic_consts[0][5,5]/GPa
       results_dict.update({'c11' : c11, 'c22' : c22, 'c33' : c33, 'c12': c12, 'c13' : c13, 'c23' : c23,
                            'c44' : c44, 'c55' : c55, 'c66' : c66})
   elif lattice_type == 'tetragonal':
       elastic_consts = matscipy.elasticity.fit_elastic_constants(bulk, symmetry='tetragonal_high', optimizer=opt, logfile=sys.stdout)
       c11 = elastic_consts[0][0,0]/GPa
       c33 = elastic_consts[0][2,2]/GPa
       c12 = elastic_consts[0][0,1]/GPa
       c13 = elastic_consts[0][0,2]/GPa
       c44 = elastic_consts[0][3,3]/GPa
       c66 = elastic_consts[0][5,5]/GPa
       results_dict.update({'c11' : c11, 'c33' : c33, 'c12': c12, 'c13' : c13, 'c44' : c44, 'c66' : c66,
                            'B' : VRH_B(c11, c33, c12, c13, c44, c66)})
   elif lattice_type == 'hexagonal':
       # Need to check if hexagonal structures are truly trigonal_high
       # symmetry=triginal_high not hexagonal until matscipy is debugged
       elastic_consts = matscipy.elasticity.fit_elastic_constants(bulk, symmetry='trigonal_high', optimizer=opt, logfile=sys.stdout)
       c11 = elastic_consts[0][0,0]/GPa
       c33 = elastic_consts[0][2,2]/GPa
       c12 = elastic_consts[0][0,1]/GPa
       c13 = elastic_consts[0][0,2]/GPa
       c44 = elastic_consts[0][3,3]/GPa
       c14 = elastic_consts[0][0,3]/GPa
       c15 = elastic_consts[0][0,4]/GPa
       c25 = elastic_consts[0][1,4]/GPa
       c66 = elastic_consts[0][5,5]/GPa
       results_dict.update({'c11' : c11, 'c33' : c33, 'c12': c12, 'c13' : c13, 'c44' : c44, 'c14' : c14,
                            'c15' : c15, 'c25' : c25, 'c66' : c66, 'B' : HTT_B(c11, c33, c12, c13)})
   elif lattice_type == 'trigonal':
       elastic_consts = matscipy.elasticity.fit_elastic_constants(bulk, symmetry='trigonal_high', optimizer=opt, logfile=sys.stdout)
       c11 = elastic_consts[0][0,0]/GPa
       c33 = elastic_consts[0][2,2]/GPa
       c12 = elastic_consts[0][0,1]/GPa
       c13 = elastic_consts[0][0,2]/GPa
       c44 = elastic_consts[0][3,3]/GPa
       c14 = elastic_consts[0][0,3]/GPa
       c15 = elastic_consts[0][0,4]/GPa
       c25 = elastic_consts[0][1,4]/GPa
       c66 = elastic_consts[0][5,5]/GPa
       results_dict.update({'c11' : c11, 'c33' : c33, 'c12': c12, 'c13' : c13, 'c44' : c44, 'c14' : c14,
                            'c15' : c15, 'c25' : c25, 'c66' : c66, 'B' : HTT_B(c11, c33, c12, c13)})

   if hasattr(model, "fix_cell_dependence"):
       model.fix_cell_dependence()

   return results_dict
