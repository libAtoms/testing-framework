#!/usr/bin/env python

import json
import ase.io

####
import sys
####
def get_ref_structs(label, models, default_analysis_settings):
    ref_structs = []
    for s in default_analysis_settings.keys():
        if s.endswith("ref_struct"):
            ref_structs.append(default_analysis_settings[s])

    data = { "composition" : {}, "min_Es" : {}  }

    # read compositions

    # read minimum energies for each structure
    for model_name in models:
        data["min_Es"][model_name] = {}
        for struct_name in ref_structs:
            with open("{}model-{}-test-bulk_{}-properties.json".format(label, model_name, struct_name), "r") as f:
                d = json.load(f)
            data["min_Es"][model_name][struct_name] = min(d["E_vs_V"], key = lambda x : x[1])[1]
            at = ase.io.read("{}model-{}-test-bulk_{}-relaxed.xyz".format(label, model_name, struct_name), format="extxyz")
            for Z in set(at.get_atomic_numbers()):
                data["composition"][struct_name] = [Z, sum(at.get_atomic_numbers() == Z)]

    return data
