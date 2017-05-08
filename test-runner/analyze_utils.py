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

    models = [os.path.basename(f) for f in glob.glob(os.path.join('..', 'models',args.label,args.models_re))]
    tests = []
    for tests_re in args.tests_re.split(','):
        tests.extend( [os.path.basename(f) for f in glob.glob(os.path.join('..', 'tests',args.label,tests_re))] )

    if args.label is None:
        args.label = ""
    else:
        args.label = args.label+"-"
    return (args, models, tests, default_analysis_settings)
