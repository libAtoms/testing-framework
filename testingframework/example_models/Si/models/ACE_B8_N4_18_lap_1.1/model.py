import pyjulip
import os
import glob

model_dir = os.path.dirname(os.path.realpath(__file__))
pot_name = glob.glob(os.path.join(model_dir,"*.json"))[0]

calculator = pyjulip.ACE(pot_name)

no_checkpoint = True

name = "ACE"
