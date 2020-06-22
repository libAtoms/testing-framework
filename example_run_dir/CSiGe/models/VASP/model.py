# Model using VASP

import sys, os
import numpy as np

from ase.calculators.vasp import Vasp2

mydir=os.path.abspath(os.path.dirname(__file__))

# assume PAWs are saved in mydir/potpaw_PBE/<element>/POTCAR
# note that different XC functionals may lead to potpaw_PBE/ replaced with another directory
# see Vasp2 calculator for details
os.environ["VASP_PP_PATH"]=mydir

# VASP KSPACING units
kspacing_value=0.15
# VASP ENCUT units
encut_value=238.0

# System-specific things below are XC functional, elec. temperature, SCF algo details chosen for 
# this material, some runtime-parallelization things like NCORE and KPAR.
# They should probably be rearranged more sensibly, some possibly externally controlled (e.g. NCORE
# maybe via env variable set by job that runs test)
#####################################################################################################

# wipe restart files (when going to new atoms object) so VASP isn't confused
# see vasp2_set_atoms() below
def wipe_restart(directory):
    try:
        os.unlink(os.path.join(directory,"WAVECAR"))
        sys.stderr.write("wiped {}\n".format(os.path.join(directory,"WAVECAR")))
    except FileNotFoundError:
        pass
    try:
        os.unlink(os.path.join(directory,"CHGCAR"))
        sys.stderr.write("wiped {}\n".format(os.path.join(directory,"CHGCAR")))
    except FileNotFoundError:
        pass

# default calculator to kspacing_value, gamma-centered
default_calculator = Vasp2(encut=encut_value,kspacing=kspacing_value,kgamma=True,
    xc='PBE', prec='Accurate',ismear=0,sigma=0.05,ediff=1.0e-8,nelm=150,
    ignore_constraints=True, ncore=16,lscalapack=False,lplane=False,
    addgrid=True,lreal=False,
    algo='normal', amix=0.03, bmix=0.01,
    isym=0,isif=3, kpar=8)

calculator = default_calculator
# wipe KPOINTS in case running in old dir and first run might be KSPACING based
try:
    os.unlink(os.path.join(calculator.directory,"KPOINTS"))
except FileNotFoundError:
    pass
wipe_restart(calculator.directory)

# monkey patch Vasp2.set_atoms to wipe restart files if major changes found
def vasp2_set_atoms(self, atoms):
    sys.stderr.write("monkey-patched set_atoms checking\n")
    if hasattr(self, "atoms") and atoms != self.atoms:
        sys.stderr.write("monkey-patched set_atoms trying to wipe\n")
        self.results = {}
        wipe_restart(self.directory)
    self.atoms = atoms.copy()
Vasp2.set_atoms = vasp2_set_atoms

# explicitly set number of k-points from now on using kspacing_value, still gamma-centered
# driver will write KPOINTS file
# at=None undoes fixed k-points
def fix_cell_dependence(at=None):
    global calculator
    if at is None:
        calculator = default_calculator
        try:
            os.unlink(os.path.join(calculator.directory,"KPOINTS"))
        except:
            pass
        print("fix_cell_dependence() going back to default")
    else:
        bz_cell = at.get_reciprocal_cell()
        n_kpts = np.floor(np.linalg.norm(bz_cell,axis=1)*2.0*np.pi/kspacing_value).astype(int)+1
        print("fix_cell_dependence() got n_kpts", n_kpts,"from recip cell",bz_cell)
        calculator = Vasp2(encut=encut_value,kpts=n_kpts,gamma=True,
            xc='PBE', prec='Accurate',ismear=0,sigma=0.05,ediff=1.0e-8,nelm=150,
            ignore_constraints=True, ncore=16,lscalapack=False,lplane=False,
            addgrid=True,lreal=False,
            algo='normal', amix=0.03, bmix=0.01,
            isym=0,isif=3, kpar=8)

# want checkpointing, but broken now because of force_consistent
# no_checkpoint = False
no_checkpoint = True
name = 'VASP'
