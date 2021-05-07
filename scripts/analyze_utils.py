#!/usr/bin/env python

import argparse
import glob
import itertools
import json
import math
import os
import sys

import ase.io

debug = False


# greatest common divisor
def gcd(l):
    divisor = l[0]
    for i in range(1, len(l)):
        divisor = math.gcd(divisor, l[i])
    return divisor


# formula unit from an array of atomic numbers
def formula_unit(atomic_numbers):
    composition = []
    for Z in sorted(set(atomic_numbers)):
        composition.append([Z, sum(atomic_numbers == Z)])
    divisor = gcd([x[1] for x in composition])
    return [[x[0], x[1] / divisor] for x in composition]


def get_matching_from_test_sets(test_sets, filename, unique=False):
    if isinstance(test_sets, str):
        test_sets = [test_sets]

    # list of globs
    files = []
    for ts in test_sets:
        files.extend(glob.glob(f'{ts}-{filename}'))

    if unique:
        if len(files) > 1:
            raise RuntimeError(f'unique but more than one file for {test_sets}, {filename} = {files}')
        elif len(files) == 0:
            return None
        else:
            return files[0]
    else:
        return files

# read a reference bulk structure, return lowest E and corresponding V, and composition
def read_ref_bulk_model_struct(test_set, model_name, struct_name):
    # min E from properties

    with open(get_matching_from_test_sets(test_set, f'model-{model_name}-test-{struct_name}-properties.json', unique=True) , "r") as f:
        d = json.load(f)
    min_EV = min(d["E_vs_V"], key = lambda x : x[1])[1]

    # composition from relaxed struct
    at = ase.io.read(get_matching_from_test_sets(test_set, f'model-{model_name}-test-{struct_name}-relaxed.xyz', unique=True), format="extxyz")
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
        print("get_multicomponent_constraints", multicomponent_constraints_structs, isinstance(multicomponent_constraints_structs, str))

    for model_name in models:
        if debug:
            print("get_multicomponent_constraints model_name", model_name)
        energy_data[model_name] = {}
        if isinstance(multicomponent_constraints_structs, str):
            # print("globbing multicomponent_constraints_struct")
            struct_name_list = get_matching_from_test_sets(test_set, f'model-{model_name}-test-{multicomponent_constraints_structs}-properties.json')
            struct_name_list = [ re.sub(f'.*-model-{model_name}-test-', '', x).replace('-properties.json','') for x in struct_name_list ]
        else:
            # print("not globbing multicomponent_constraints_struct")
            struct_name_list = multicomponent_constraints_structs
        # print("get_multicomponent_constraints struct_name_list",struct_name_list)

        for struct_name in struct_name_list:
            if debug:
                print("get_multicomponent_constraints struct_name",test_set, model_name, struct_name)
            try:
                (min_EV, composition) = read_ref_bulk_model_struct(test_set, model_name, struct_name)
                energy_data[model_name][struct_name] = min_EV
                if not "struct_name" in composition_data:
                    composition_data[struct_name] = composition
                if debug:
                    print("min_EV, composition ",min_EV, composition)
            except:
                pass

    return (composition_data, energy_data)

# start an analysis by parsing arguments, reading defaults
def analyze_start(default_tests=['*']):
    parser = argparse.ArgumentParser(description='Analyze bulk lattices')
    parser.add_argument('--models', '-m', action='store', nargs='+', type=str, help='models to include', default=['*'])
    if not isinstance(default_tests,list):
        default_tests=[default_tests]
    parser.add_argument('--tests', '-t', action='store', nargs='+', type=str, help='tests to include (globs)', default=default_tests)
    # duplicates of stuff in scripts/run-all.py
    parser.add_argument('--test_set', '-s', nargs='+', action='store', help='label for tests directories', required=True)
    parser.add_argument('--models_path', '-P', action='store', help='path to models directory', default=os.path.join(os.getcwd(),'../models'))
    parser.add_argument('--REF_E_offset', type=float, action='store', help='offset of reference model energy/atom', default=0.0)

    try:
        with open("default_analysis_opts.json","r") as f:
            try:
                default_analysis_opts = json.load(f)
            except Exception as e:
                sys.stderr.write("ERROR json.load on file default_analysis_opts.json\n")
                raise e
    except:
        default_analysis_opts=[]

    args = parser.parse_args(sys.argv[1:] + default_analysis_opts)

    # from full list of models/tests
    models = [ os.path.split(f)[1] for f in list(itertools.chain.from_iterable([ glob.glob(os.path.join(args.models_path, d)) for d in args.models ])) ]
    my_path=os.path.dirname(os.path.realpath(__file__))
    tests = []
    for test_set in args.test_set:
        tests_path = os.path.join(my_path,"..","tests",test_set)
        tests.extend([ os.path.split(f)[1] for f in list(itertools.chain.from_iterable([ glob.glob(os.path.join(tests_path, d)) for d in args.tests ])) ])

    with open("default_analysis_settings.json") as default_analysis_file:
        try:
            default_analysis_settings = json.load(default_analysis_file)
        except Exception as e:
            sys.stderr.write("ERROR json.load on file default_analysis_settings.json\n")
            raise e

    return (args, models, tests, default_analysis_settings)

# read the properties for a set of models and tests with a given system
def read_properties(models, tests, test_set):
    data = {}
    for model_name in models:
        print("reading data for model {}".format(model_name))
        data[model_name] = {}
        cur_model_data = {}
        for test_name in tests:
            print("   reading data for test {}".format(test_name))
            try:
                prop_filename = get_matching_from_test_sets(test_set, f'model-{model_name}-test-{test_name}-properties.json', unique=True)
            except:
                print("Failed to get unique prop_filename for ", test_set, f'model-{model_name}-test-{test_name}-properties.json')
                continue

            try:
                with open(prop_filename, "r") as model_data_file:
                    cur_model_data[test_name] = json.load(model_data_file)
            except:
                print("No properties file '{}'".format(prop_filename))
                continue

        data[model_name] = cur_model_data

    return data

