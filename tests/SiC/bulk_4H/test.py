import lattice
import os.path

properties = lattice.do_lattice(os.path.abspath(os.path.dirname(__file__)), 'hexagonal')
