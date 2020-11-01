import os.path

import lattice

properties = lattice.do_lattice(os.path.abspath(os.path.dirname(__file__)),
                                'hexagonal')  # , method='sd2') - for Brenner # To be fixed soon
