# This script defines a test case which computes one or more physical
# properties with a given model
#
# INPUTS: 
#   model.calculator -- an ase.calculator.Calculator instance
#     this script can assume the calculator is checkpointed.
#
# OUTPUTS:
#   properties -- dictionary of key/value pairs corresponding
#     to physical quantities computed by this test

# standard ASE structure generation routines
import ase, os

import lattice_cubic

# the current model
import model 

test_dir = os.path.abspath(os.path.dirname(__file__))

# set up the a
bulk = ase.io.read(test_dir+"/bulk.xyz", format="extxyz")

(c11, c12, c44, E_vs_V) = lattice_cubic.do_lattice(bulk, __file__)
B = (c11+2.0*c12)/3.0

# dictionary of computed properties - this is output of this test, to
#   be compared with other models
properties = {'c11' : c11, 'c12' : c12, 'c44' : c44, 'bulk_modulus' : B, 'E_vs_V' : E_vs_V }
