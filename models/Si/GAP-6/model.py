# Model for GAP-6 (Si GAP to be published in PRX)

from quippy import Potential
import os.path

calculator = Potential(param_filename=os.path.join(os.path.abspath(os.path.dirname(__file__)),'gp_iter6_sparse9k.xml'))

no_checkpoint = True
name = 'GAP-6'
