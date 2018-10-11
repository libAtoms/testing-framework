# Model for VASP for Si

from ase.calculators.vasp import Vasp
from ase.dft.kpoints import monkhorst_pack
import os.path
import numpy as np

# define calculator symbol

mydir=os.path.abspath(os.path.dirname(__file__))
os.environ["VASP_PP_PATH"]=mydir
kspacing_value=0.250

calculator = Vasp(encut=400.0,kpts=1,
    xc='pbesol',prec='Accurate',ismear=0,sigma=0.05,ediff=1.0e-7,nelm=150,
    ncore=16,lscalapack=False,lplane=False,
    addgrid=True,lreal=False,
    algo='normal',amix=0.1,
    isym=0,isif=3)
# calculator.debug = 'irreducible'

def new_cell(at):
    global calculator
    bz_cell = at.get_reciprocal_cell()
    n_kpts = np.floor(np.linalg.norm(bz_cell,axis=1)*2.0*np.pi/kspacing_value).astype(int)+1
    print "new_cell got n_kpts", n_kpts
    calculator = Vasp(encut=400.0,kpts=n_kpts,gamma=True,
        xc='pbesol',prec='Accurate',ismear=0,sigma=0.05,ediff=1.0e-7,nelm=150,
        ncore=16,lscalapack=False,lplane=False,
        addgrid=True,lreal=False,
        algo='normal',amix=0.1,
        isym=0,isif=3)
    # calculator.debug = 'irreducible'

# want checkpointing, but broken now because of force_consistent
# no_checkpoint = False
no_checkpoint = True
name = 'VASP'
