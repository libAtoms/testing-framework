# Model for GAP-6 (Si GAP to be published in PRX)

from quippy import Potential
import os.path

calculator = Potential(param_filename=os.path.join(os.path.abspath(os.path.dirname(__file__)),'GAP.iter_6.xml'))

no_checkpoint = True
name = 'gap-rss-test-1'
