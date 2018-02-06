import numpy as np
from utilities import relax_config, run_root
import ase.io, sys, os.path
from ase.optimize.precon import Exp, PreconLBFGS
import matscipy.elasticity
from ase.units import GPa

# rotate to align (111) with z
def align_111_with_z(bulk):
    lat = bulk.get_cell()

    old_z = np.array(lat[0,:] + lat[1,:] + lat[2,:]); old_z /= np.linalg.norm(old_z)
    old_a1 = lat[0,:]
    a1_ctr_angle = np.arccos(np.dot(old_z, old_a1/np.linalg.norm(old_a1)))
    new_a1 = np.array([1.0,0.0,0.0])*np.sin(a1_ctr_angle) + np.array([0.0,0.0,1.0])*np.cos(a1_ctr_angle)
    new_other_dir = np.cross(old_z, old_a1); new_other_dir /= np.linalg.norm(new_other_dir)

    old_lat = np.zeros( (3,3) )
    new_lat = np.zeros( (3,3) )
    old_lat[:,0] = old_z # z
    old_lat[:,1] = old_a1/np.linalg.norm(old_a1) # a1
    old_lat[:,2] = new_other_dir # other dir
    new_lat[:,0] = [0.0, 0.0, 1.0]
    new_lat[:,1] = new_a1
    new_lat[:,2] = [0.0, 1.0, 0.0]

    # new_lat = m.old_lat
    # m = new_lat . old_lat_inv
    m = np.dot(new_lat, np.linalg.inv(old_lat))

    new_cell = np.dot(m, lat.T).T

    bulk.set_cell(new_cell, scale_atoms=True)

# Mater Trans v 53 p 1247 (2012), Eq. 4-10
def VRH_B(c11, c33, c12, c13, c44, c66):
    Bv = (1.0/9.0)*(2.0*(c11+c12) + c33 + 4.0*c13)

    M = c11 + c12 + 2.0*c33 - 4.0*c13
    Csq = (c11+c12)*c33 - 2.0*c13**2
    Br = Csq/M

    return (Bv+Br)/2.0

def calc_E_vs_V(bulk, vol_range=0.25, n_steps=10, tol=1e-2, method='lbfgs'): # hack tol to deal with Te C2/m
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
          scaled_bulk = relax_config(scaled_bulk, relax_pos=True, relax_cell=True, tol=tol, traj_file=None, constant_volume=True, method=method,
              refine_symmetry_tol=1.0e-4, keep_symmetry=True, config_label="E_vs_V_%02d" % i, from_base_model=True, save_config=True)
      except:
          print "failed"
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
          scaled_bulk = relax_config(scaled_bulk, relax_pos=True, relax_cell=True, tol=tol, traj_file=None, constant_volume=True, method=method,
              refine_symmetry_tol=1.0e-4, keep_symmetry=True, config_label="E_vs_V_%02d" % i, from_base_model=True, save_config=True)
      except:
          print "failed"
          break
      ase.io.write(sys.stdout, scaled_bulk, format='extxyz')
      E_vs_V.append( (scaled_bulk.get_volume()/len(scaled_bulk), scaled_bulk.get_potential_energy()/len(bulk)) )

   return E_vs_V


def do_lattice(test_dir, lattice_type, vol_range=0.25):

   bulk = ase.io.read(test_dir+"/bulk.xyz", format="extxyz")

   tol = 1e-3 # max force tol for relaxation

   print "relax bulk"
   # relax the initial unit cell and atomic positions
   bulk = relax_config(bulk, relax_pos=True, relax_cell=True, tol=tol, traj_file=None, method='cg_n', 
     refine_symmetry_tol=1.0e-2, keep_symmetry=True, config_label="bulk", from_base_model=True, save_config=True)

   print "final relaxed bulk"
   ase.io.write(sys.stdout, bulk, format='extxyz')
   ase.io.write(os.path.join("..",run_root+"-relaxed.xyz"),  bulk, format='extxyz')

   print "calculating E vs. V"
   E_vs_V = calc_E_vs_V(bulk, vol_range=vol_range)

   print "calculating elastic constants"
   precon = Exp(3.0)
   opt = lambda atoms, **kwargs: PreconLBFGS(atoms, precon=precon, **kwargs)
   if lattice_type == 'cubic':
       elastic_consts = matscipy.elasticity.fit_elastic_constants(bulk, symmetry='cubic', optimizer=opt, logfile=sys.stdout)
       c11 = elastic_consts[0][0,0]/GPa
       c12 = elastic_consts[0][0,1]/GPa
       c44 = elastic_consts[0][3,3]/GPa
       return ({'c11' : c11, 'c12': c12, 'c44' : c44, 'B' : (c11+2.0*c12)/3.0, 'E_vs_V' : E_vs_V})
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
       return ({'c11' : c11, 'c22' : c22, 'c33' : c33, 'c12': c12, 'c13' : c13, 'c23' : c23,
           'c44' : c44, 'c55' : c55, 'c66' : c66, 'E_vs_V' : E_vs_V})
   elif lattice_type == 'tetragonal':
       elastic_consts = matscipy.elasticity.fit_elastic_constants(bulk, symmetry='tetragonal_high', optimizer=opt, logfile=sys.stdout)
       c11 = elastic_consts[0][0,0]/GPa
       c33 = elastic_consts[0][2,2]/GPa
       c12 = elastic_consts[0][0,1]/GPa
       c13 = elastic_consts[0][0,2]/GPa
       c44 = elastic_consts[0][3,3]/GPa
       c66 = elastic_consts[0][5,5]/GPa
       return ({'c11' : c11, 'c33' : c33, 'c12': c12, 'c13' : c13, 'c44' : c44, 'c66' : c66,
        'B' : VRH_B(c11, c33, c12, c13, c44, c66), 'E_vs_V' : E_vs_V})
