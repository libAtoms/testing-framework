import numpy as np
from utilities import relax_config, run_root
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

def calc_E_vs_V(bulk, vol_range=0.25, n_steps=10, tol=1e-2, method='lbfgs'): # hack tol to deal with Te C2/m
   import model

   V0 = bulk.get_volume()
   dV = bulk.get_volume()*vol_range/n_steps
   E_vs_V=[]

   scaled_bulk = bulk.copy()
   for i in range(0, -n_steps-1, -1):
      V_cur = scaled_bulk.get_volume()
      scaled_bulk.set_cell(scaled_bulk.get_cell()*((V0+i*dV)/V_cur)**(1.0/3.0), scale_atoms=True)
      ase.io.write(sys.stdout, scaled_bulk, format='extxyz')
      print "trying to relax i",i
      try:
          if hasattr(model, "fix_cell_dependence"):
               model.fix_cell_dependence(scaled_bulk)
          scaled_bulk = relax_config(scaled_bulk, relax_pos=True, relax_cell=True, tol=tol, traj_file=None, constant_volume=True, method=method,
              refine_symmetry_tol=1.0e-4, keep_symmetry=True, config_label="E_vs_V_%02d" % i, from_base_model=True, save_config=True)
      except Exception, e:
          print "WARNING: failed config in calc_E_vs_V", str(e)
          break
      ase.io.write(sys.stdout, scaled_bulk, format='extxyz')
      E_vs_V.insert(0, (scaled_bulk.get_volume()/len(scaled_bulk), scaled_bulk.get_potential_energy()/len(bulk)) )

   scaled_bulk = bulk.copy()
   for i in range(1,n_steps+1):
      V_cur = scaled_bulk.get_volume()
      scaled_bulk.set_cell(scaled_bulk.get_cell()*((V0+i*dV)/V_cur)**(1.0/3.0), scale_atoms=True)
      ase.io.write(sys.stdout, scaled_bulk, format='extxyz')
      print "trying to relax i",i
      try:
          if hasattr(model, "fix_cell_dependence"):
               model.fix_cell_dependence(scaled_bulk)
          scaled_bulk = relax_config(scaled_bulk, relax_pos=True, relax_cell=True, tol=tol, traj_file=None, constant_volume=True, method=method,
              refine_symmetry_tol=1.0e-4, keep_symmetry=True, config_label="E_vs_V_%02d" % i, from_base_model=True, save_config=True)
      except Exception, e:
          print "failed", str(e)
          break
      ase.io.write(sys.stdout, scaled_bulk, format='extxyz')
      E_vs_V.append( (scaled_bulk.get_volume()/len(scaled_bulk), scaled_bulk.get_potential_energy()/len(bulk)) )

   if hasattr(model, "fix_cell_dependence"):
       model.fix_cell_dependence()

   return E_vs_V


def do_lattice(test_dir, lattice_type, vol_range=0.25, method='lbfgs'):

   import model
   bulk = ase.io.read(test_dir+"/bulk.xyz", format="extxyz")

   results_dict = {}

   tol = 1e-2 # max force tol for relaxation

   print "relax bulk"
   # relax the initial unit cell and atomic positions
   if hasattr(model, "fix_cell_dependence"):
       model.fix_cell_dependence(bulk)
   bulk = relax_config(bulk, relax_pos=True, relax_cell=True, tol=tol, traj_file="lattice_bulk_traj.xyz", method=method, 
                       refine_symmetry_tol=1.0e-2, keep_symmetry=True, config_label="bulk", from_base_model=True, save_config=True)
   if hasattr(model, "fix_cell_dependence"):
       model.fix_cell_dependence()

   print "final relaxed bulk"
   ase.io.write(sys.stdout, bulk, format='extxyz')
   ase.io.write(os.path.join("..",run_root+"-relaxed.xyz"),  bulk, format='extxyz')

   print "calculating E vs. V"
   E_vs_V = calc_E_vs_V(bulk, vol_range=vol_range)
   results_dict.update({ 'E_vs_V' : E_vs_V })

   print "calculating elastic constants"

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
