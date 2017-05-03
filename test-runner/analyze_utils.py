#!/usr/bin/env python

import json
import ase.io

def get_element_ref_structs(label, models, element_ref_struct):
    data = {}

    for symb in element_ref_struct:
        struct_name = element_ref_struct[symb]
        data[symb] = { "min_Es" : {}, "composition" : {}, "ref_struct" : struct_name }
        for model_name in models:
            data[symb]["min_Es"][model_name] = {}
            with open("{}model-{}-test-bulk_{}-properties.json".format(label, model_name, struct_name), "r") as f:
                d = json.load(f)
            data[symb]["min_Es"][model_name] = min(d["E_vs_V"], key = lambda x : x[1])[1]
            at = ase.io.read("{}model-{}-test-bulk_{}-relaxed.xyz".format(label, model_name, struct_name), format="extxyz")
            for Z in set(at.get_atomic_numbers()):
                data[symb]["composition"] = [Z, sum(at.get_atomic_numbers() == Z)]

    return data
