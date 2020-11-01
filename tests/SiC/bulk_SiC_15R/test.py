import os.path

import lattice

properties = lattice.do_lattice(os.path.abspath(os.path.dirname(__file__)), 'trigonal')
