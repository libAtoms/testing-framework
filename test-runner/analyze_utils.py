#!/usr/bin/env python

import json
import ase.io
import fractions
import argparse, os, glob

def gcd(l):
    divisor = l[0]
    for i in range(1, len(l)):
        divisor = fractions.gcd(divisor, l[i])
    return divisor

def formula_unit(atomic_numbers):
    composition = []
    for Z in set(atomic_numbers):
        composition.append([Z, sum(atomic_numbers == Z)])
    divisor = gcd([ x[1] for x in composition ])
    return [ [x[0], x[1]/divisor] for x in composition ]

def read_ref_bulk_model_struct(label, model_name, struct_name=None, struct_file=None):
    # min E from properties
    with open("{}model-{}-test-{}-properties.json".format(label, model_name, struct_name), "r") as f:
        d = json.load(f)
    min_EV = min(d["E_vs_V"], key = lambda x : x[1])[1]

    # composition from relaxed struct
    at = ase.io.read("{}model-{}-test-{}-relaxed.xyz".format(label, model_name, struct_name), format="extxyz")
    composition = formula_unit(at.get_atomic_numbers())

    return (min_EV, composition)

def get_element_ref_structs(label, models, element_ref_struct):
    data = {}

    for symb in element_ref_struct:
        struct_name = element_ref_struct[symb]
        data[symb] = { "min_Es" : {}, "ref_struct" : struct_name }
        for model_name in models:
            try:
                (min_EV, composition) = read_ref_bulk_model_struct(label, model_name, struct_name)
            except:
                continue
            data[symb]["min_Es"][model_name] = min_EV
            if not "composition" in data[symb]: 
                data[symb]["composition"] = composition
    return data

def get_multicomponent_constraints(label, models, multicomponent_constraints_structs):
    composition_data = {}
    energy_data = {}

    for model_name in models:
        # print "get_multicomponent_constraints", struct_name
        energy_data[model_name] = {}
        if istype(multicomponent_constraints_structs, str):
            struct_name_list = glob.glob("{}model-{}-test-{}-properties.json".format(label, model_name, multicomponent_constraints_structs))
            struct_name_list = [ x.replace("{}model-{}-test-".format(label.model_name),"").replace("-properties.json","") for x in struct_name_list ]
        else:
            struct_name_list = multicomponent_constraints_structs

        print "get_multicomponent_constraints", struct_name_list
        sys.exit(1)

        for struct_name in struct_name_list:
            (min_EV, composition) = read_ref_bulk_model_struct(label, model_name, struct_name)
            energy_data[model_name][struct_name] = min_EV
            if not "struct_name" in composition_data: 
                composition_data[struct_name] = composition

    return (composition_data, energy_data)

def analyze_start(default_tests_re):
    parser = argparse.ArgumentParser(description='Analyze bulk lattices')
    parser.add_argument('--models_re', '-m', action='store', type=str, help='models to include', default='*')
    parser.add_argument('--tests_re', '-t', action='store', type=str, help='tests to include', default=default_tests_re)
    parser.add_argument('--label', '-l', action='store', help='optional label for models/tests directories', default='')
    args = parser.parse_args()

    try:
        with open("DEFAULTS_LABEL","r") as f:
            defaults_label = f.readline().strip()+"_"
    except:
        defaults_label = ""

    with open("{}default_analysis_settings.json".format(defaults_label)) as default_analysis_file:
        default_analysis_settings = json.load(default_analysis_file)

    with open("{}default_run_opts.json".format(defaults_label)) as default_run_file:
        default_run_opts = json.load(default_run_file)
        if 'label' in default_run_opts and not args.label:
            args.label = default_run_opts['label']

    models = []
    for models_re in args.models_re.split(','):
        models.extend( [os.path.basename(f) for f in glob.glob(os.path.join('..', 'models',args.label,models_re))] )
    tests = []
    for tests_re in args.tests_re.split(','):
        tests.extend( [os.path.basename(f) for f in glob.glob(os.path.join('..', 'tests',args.label,tests_re))] )

    if args.label is None:
        args.label = ""
    else:
        args.label = args.label+"-"
    return (args, models, tests, default_analysis_settings)

def read_properties(models, tests, label):
    data = {}
    for model_name in models:
        print "reading data for model {}".format(model_name)
        data[model_name] = {}
        cur_model_data = {}
        for test_name in tests:
            print "   reading data for test {}".format(test_name)

            prop_filename ="{}model-{}-test-{}-properties.json".format(label, model_name, test_name)
            try:
                with open(prop_filename, "r") as model_data_file:
                    cur_model_data[test_name] = json.load(model_data_file)
            except:
                print "No properties file '{}'".format(prop_filename)
                continue

        data[model_name] = cur_model_data

    return data

