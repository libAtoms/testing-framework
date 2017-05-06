#!/usr/bin/env python

import json
import ase.io
import fractions

def gcd(l):
    divisor = l[0]
    for i in range(1, len(l)):
        divisor = fractions.gcd(divisor, l[i])
    return divisor

def read_ref_bulk_model_struct(label, model_name, struct_name):
    # min E from properties
    with open("{}model-{}-test-{}-properties.json".format(label, model_name, struct_name), "r") as f:
        d = json.load(f)
    min_EV = min(d["E_vs_V"], key = lambda x : x[1])[1]

    # composition from relaxed struct
    composition = []
    at = ase.io.read("{}model-{}-test-{}-relaxed.xyz".format(label, model_name, struct_name), format="extxyz")
    for Z in set(at.get_atomic_numbers()):
        composition.append([Z, sum(at.get_atomic_numbers() == Z)])

    divisor = gcd([ x[1] for x in composition ])
    composition = [ [x[0], x[1]/divisor] for x in composition ]

    return (min_EV, composition)

def get_element_ref_structs(label, models, element_ref_struct):
    data = {}

    for symb in element_ref_struct:
        struct_name = element_ref_struct[symb]
        data[symb] = { "min_Es" : {}, "ref_struct" : struct_name }
        for model_name in models:
            (min_EV, composition) = read_ref_bulk_model_struct(label, model_name, struct_name)
            data[symb]["min_Es"][model_name] = min_EV
            if not "composition" in data[symb]: 
                data[symb]["composition"] = composition
    return data

def get_multicomponent_constraints(label, models, multicomponent_constraints):
    composition_data = {}
    energy_data = {}

    for model_name in models:
        # print "get_multicomponent_constraints", struct_name
        energy_data[model_name] = {}
        for struct_name in multicomponent_constraints:
            (min_EV, composition) = read_ref_bulk_model_struct(label, model_name, struct_name)
            energy_data[model_name][struct_name] = min_EV
            if not "struct_name" in composition_data: 
                composition_data[struct_name] = composition

    return (composition_data, energy_data)
