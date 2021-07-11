# original ANI model

import torch
from torchani.neurochem import load_model
from torchani.ase import Calculator
from torchani import AEVComputer
from torchani.neurochem import load_sae
from torchani.nn import Sequential
import os

# A module defining a module needs to define only one variable,
# named `calculator`, which should be an instance of the ase.calculator.Calculator,
# a subclass of this, or a compatible class implementing the calculator interface.

device = torch.device('cuda' if torch.cuda.is_available() else 'cpu')
model_dir = os.path.dirname(os.path.realpath(__file__))

species_order = ['H', 'C', 'N', 'O']
const_file = os.path.join(model_dir, '../ANI_common/rHCNO-5.2R_16-3.5A_a4-8.params')
consts = torchani.neurochem.Constants(const_file)
aev_computer = AEVComputer(**consts)
energy_shifter = load_sae(os.path.join(model_dir, '../ANI_common/sae_linfit.dat'))


nn = load_model(consts.species, os.path.join(model_dir, '../ANI_common'))
model_trained = Sequential(aev_computer, nn, energy_shifter).to(device)

calculator = Calculator(species=species_order, model=model_trained)

name = 'ANI_orig'
