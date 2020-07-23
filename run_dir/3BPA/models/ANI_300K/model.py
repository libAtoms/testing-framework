# ANI model fitted on 300K configs 

from torchani.neurochem import load_model
from torchani.ase import Calculator
from torchani import AEVComputer
from torchani.neurochem import load_sae
from torchani.nn import Sequential

# A module defining a module needs to define only one variable,
# named `calculator`, which should be an instance of the ase.calculator.Calculator,
# a subclass of this, or a compatible class implementing the calculator interface.

species_order = ['H', 'C', 'N', 'O']
const_file = '../ANI_common/rHCNO-5.2R_16-3.5A_a4-8.params')
consts = torchani.neurochem.Constants(const_file)
aev_computer = AEVComputer(**consts)
energy_shifter = load_sae('../ANI_common/sae_linfit.dat')



nn = load_model(consts.species, '../ANI_common')
nn.load_state_dict(torch.load('model300EF_pre_2_best.pt'))
model_trained = Sequential(aev_computer, nn, energy_shifter).to(device)

calculator = Calculator(species=species_order, model=model_trained)

name = 'ANI_300K'
