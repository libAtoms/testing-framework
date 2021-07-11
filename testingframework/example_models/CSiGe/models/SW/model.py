# Model for Stillinger-Weber with original parameters for Si (Z=14)

from quippy.potential import Potential
from utilities import path_of_file


# A module defining a module needs to define only one variable,
# named `calculator`, which should be an instance of the ase.calculator.Calculator,
# a subclass of this, or a compatible class implementing the calculator interface.

calculator = Potential('IP SW', param_filename=(path_of_file(__file__)+"/ip.parms.SW.xml"))
no_checkpoint = True

name = 'SW'
