from atomistica import Tersoff, Tersoff_PRB_39_5566_Si_C

# A module defining a module needs to define only one variable,
# named `calculator`, which should be an instance of the ase.calculator.Calculator,
# a subclass of this, or a compatible class implementing the calculator interface.

calculator = Tersoff(**Tersoff_PRB_39_5566_Si_C)

no_checkpoint = True
name = 'Tersoff'
