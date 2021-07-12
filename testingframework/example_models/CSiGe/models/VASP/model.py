# Model using VASP

import sys, os
import numpy as np

from ase.calculators.vasp import Vasp2

mydir = os.path.abspath(os.path.dirname(__file__))

# assume PAWs are saved in mydir/potpaw_PBE/<element>/POTCAR
# note that different XC functionals may lead to potpaw_PBE/ replaced with another directory
# see Vasp2 calculator for details
os.environ["VASP_PP_PATH"] = mydir

# VASP KSPACING units
kspacing_value = 0.15

# keywords related to this model (except k-points, which are handled separately)
vasp_keywords = {
    "encut": 238.0,
    "xc": "PBE",
    "ismear": 0,
    "sigma": 0.05,
    "nelm": 150,
    "algo": "normal",
    "amix": 0.03,
    "bmix": 0.01,
}
# keywords related to accuracy
vasp_keywords.update(
    {"ediff": 1.0e-8, "prec": "Accurate", "addgrid": True, "lreal": False}
)
# keywords that should always be there
vasp_keywords.update({"ignore_constraints": True, "isif": 3, "isym": 0})
# keywords related to parallelism.  Setting ncore sensibly is important for efficiency
vasp_keywords.update(
    {"lscalapack": False, "lplane": False, "kpar": 8, "ncore": os.environ["VASP_NCORE"]}
)

# wipe restart files (when going to new atoms object) so VASP isn't confused
# see vasp2_set_atoms() below
def wipe_restart(directory):
    try:
        os.unlink(os.path.join(directory, "WAVECAR"))
        sys.stderr.write("wiped {}\n".format(os.path.join(directory, "WAVECAR")))
    except FileNotFoundError:
        pass
    try:
        os.unlink(os.path.join(directory, "CHGCAR"))
        sys.stderr.write("wiped {}\n".format(os.path.join(directory, "CHGCAR")))
    except FileNotFoundError:
        pass


# default calculator to kspacing_value, gamma-centered
default_calculator = Vasp2(kspacing=kspacing_value, kgamma=True, **vasp_keywords)
calculator = default_calculator

# wipe KPOINTS in case running in old dir and first run might be KSPACING based
try:
    os.unlink(os.path.join(calculator.directory, "KPOINTS"))
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
            os.unlink(os.path.join(calculator.directory, "KPOINTS"))
        except:
            pass
        print("fix_cell_dependence() going back to default")
    else:
        bz_cell = at.get_reciprocal_cell()
        n_kpts = (
            np.floor(
                np.linalg.norm(bz_cell, axis=1) * 2.0 * np.pi / kspacing_value
            ).astype(int)
            + 1
        )
        print("fix_cell_dependence() got n_kpts", n_kpts, "from recip cell", bz_cell)
        calculator = Vasp2(kpts=n_kpts, gamma=True, **vasp_keywords)


# want checkpointing, but broken now because of force_consistent
# no_checkpoint = False
no_checkpoint = True
name = "VASP"
