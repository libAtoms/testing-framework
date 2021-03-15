"""
Silagraphene, layered SiC structure with AA stacking of planes.

- initial structure from SI in https://doi.org/10.1021/acs.chemmater.8b03293
- relaxed with CASTEP (PBE, 500eV, 0.05 kp-spacing)

"""
import os.path

import lattice

properties = lattice.do_lattice(os.path.abspath(os.path.dirname(__file__)), 'hexagonal')
