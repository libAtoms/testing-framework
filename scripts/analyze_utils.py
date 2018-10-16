#!/usr/bin/env python

import json
import ase.io
import fractions
import itertools
import argparse, os, glob, sys

debug = False

# greatest common divisor
def gcd(l):
    divisor = l[0]
    for i in range(1, len(l)):
        divisor = fractions.gcd(divisor, l[i])
    return divisor

# formula unit from an array of atomic numbers
def formula_unit(atomic_numbers):
    composition = []
    for Z in set(atomic_numbers):
        composition.append([Z, sum(atomic_numbers == Z)])
    divisor = gcd([ x[1] for x in composition ])
    return [ [x[0], x[1]/divisor] for x in composition ]

# read a reference bulk structure, return lowest E and corresponding V, and composition
def read_ref_bulk_model_struct(test_set, model_name, struct_name):
    # min E from properties
    with open("{}-model-{}-test-{}-properties.json".format(test_set, model_name, struct_name), "r") as f:
        d = json.load(f)
    min_EV = min(d["E_vs_V"], key = lambda x : x[1])[1]

    # composition from relaxed struct
    at = ase.io.read("{}-model-{}-test-{}-relaxed.xyz".format(test_set, model_name, struct_name), format="extxyz")
    composition = formula_unit(at.get_atomic_numbers())

    return (min_EV, composition)

# get energies and volumes for every element's reference structure, for all models
def get_element_ref_structs(test_set, models, element_ref_struct):
    data = {}

    for symb in element_ref_struct:
        struct_name = element_ref_struct[symb]
        data[symb] = { "min_Es" : {}, "ref_struct" : struct_name }
        for model_name in models:
            try:
                (min_EV, composition) = read_ref_bulk_model_struct(test_set, model_name, struct_name)
            except:
                continue
            data[symb]["min_Es"][model_name] = min_EV
            if not "composition" in data[symb]: 
                data[symb]["composition"] = composition
    return data

# get the multicomponent constraints, from a list or a glob, and store their compositions and energies
def get_multicomponent_constraints(test_set, models, multicomponent_constraints_structs):
    composition_data = {}
    energy_data = {}
    if debug:
        print "get_multicomponent_constraints", multicomponent_constraints_structs, isinstance(multicomponent_constraints_structs, basestring)

    for model_name in models:
        if debug:
            print "get_multicomponent_constraints model_name", model_name
        energy_data[model_name] = {}
        if isinstance(multicomponent_constraints_structs, basestring):
            # print "globbing multicomponent_constraints_struct"
            struct_name_list = glob.glob("{}-model-{}-test-{}-properties.json".format(test_set, model_name, multicomponent_constraints_structs))
            struct_name_list = [ x.replace("{}-model-{}-test-".format(test_set, model_name),"").replace("-properties.json","") for x in struct_name_list ]
        else:
            # print "not globbing multicomponent_constraints_struct"
            struct_name_list = multicomponent_constraints_structs
        # print "get_multicomponent_constraints struct_name_list",struct_name_list

        for struct_name in struct_name_list:
            if debug:
                print "get_multicomponent_constraints struct_name",struct_name
            (min_EV, composition) = read_ref_bulk_model_struct(test_set, model_name, struct_name)
            energy_data[model_name][struct_name] = min_EV
            if not "struct_name" in composition_data: 
                composition_data[struct_name] = composition
            if debug:
                print "min_EV, composition ",min_EV, composition

    return (composition_data, energy_data)

# start an analysis by parsing arguments, reading defaults
def analyze_start(default_tests=['*']):
    parser = argparse.ArgumentParser(description='Analyze bulk lattices')
    parser.add_argument('--models', '-m', action='store', nargs='+', type=str, help='models to include', default=['*'])
    if not isinstance(default_tests,list):
        default_tests=[default_tests]
    parser.add_argument('--tests', '-t', action='store', nargs='+', type=str, help='tests to include (globs)', default=default_tests)
    # duplicates of stuff in scripts/run-all.py
    parser.add_argument('--test_set', '-s', action='store', help='label for tests directories', required=True)
    parser.add_argument('--models_path', '-P', action='store', help='path to models directory', default=os.path.join(os.getcwd(),'../models'))

    try:
        with open("default_analysis_opts.json","r") as f:
            default_analysis_opts = json.load(f)
    except:
        default_analysis_opts=[]

    args = parser.parse_args(sys.argv[1:] + default_analysis_opts)

    # from full list of models/tests
    models = [ os.path.split(f)[1] for f in list(itertools.chain.from_iterable([ glob.glob(os.path.join(args.models_path, d)) for d in args.models ])) ]
    my_path=os.path.dirname(os.path.realpath(__file__))
    tests_path = os.path.join(my_path,"..","tests",args.test_set)
    tests = [ os.path.split(f)[1] for f in list(itertools.chain.from_iterable([ glob.glob(os.path.join(tests_path, d)) for d in args.tests ])) ]

    with open("default_analysis_settings.json") as default_analysis_file:
        default_analysis_settings = json.load(default_analysis_file)

    return (args, models, tests, default_analysis_settings)

# read the properties for a set of models and tests with a given system
def read_properties(models, tests, test_set):
    data = {}
    for model_name in models:
        print "reading data for model {}".format(model_name)
        data[model_name] = {}
        cur_model_data = {}
        for test_name in tests:
            print "   reading data for test {}".format(test_name)

            prop_filename ="{}-model-{}-test-{}-properties.json".format(test_set, model_name, test_name)
            try:
                with open(prop_filename, "r") as model_data_file:
                    cur_model_data[test_name] = json.load(model_data_file)
            except:
                print "No properties file '{}'".format(prop_filename)
                continue

        data[model_name] = cur_model_data

    return data

