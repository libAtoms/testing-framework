# Model for Stillinger-Weber with original parameters for Si (Z=14)

from quippy.potential import Potential

# A module defining a module needs to define only one variable,
# named `calculator`, which should be an instance of the ase.calculator.Calculator,
# a subclass of this, or a compatible class implementing the calculator interface.

calculator = Potential('IP Tersoff', param_str="""
<Tersoff_params n_types="3" label="PRB_39">
<comment> Tersoff, Phys. Rev. B v 39, p 5566 (1989) </comment>
<per_type_data type="1" atomic_num="6" 
  A="1393.6" B="346.7" lambda="3.4879" mu="2.2119" beta="0.00000015724" 
  n="0.72751" c="38049" d="4.384" h="-0.57058" R="1.8" S="2.1" />
<per_type_data type="2" atomic_num="14"
 A="1830.8" B="471.18" lambda="2.4799" mu="1.7322" beta="0.0000011"
 n="0.78734" c="100390" d="16.217" h="-0.59825" R="2.7" S="3.0" />
<per_type_data type="3" atomic_num="32"
  A="1769" B="419.23" lambda="2.4451" mu="1.7047" beta="0.00000090166"
  n="0.75627" c="106430" d="15.652" h="-0.43884" R="2.8" S="3.1" />
<per_pair_data type1="1" type2="1" chi="1.0" />
<per_pair_data type1="2" type2="1" chi="0.9776" />
<per_pair_data type1="2" type2="2" chi="1.0" />
<per_pair_data type1="3" type2="1" chi="1.0" />
<per_pair_data type1="3" type2="2" chi="1.00061" />
<per_pair_data type1="3" type2="3" chi="1.0" />
</Tersoff_params>
""")

no_checkpoint = True

name = 'Tersoff'
