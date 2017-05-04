#!/usr/bin/env python

# import matplotlib
# matplotlib.use('PDF')
# from matplotlib.pyplot import *
import json
import argparse
import glob, os
import ase.io
from ase.data import chemical_symbols
from analyze_utils import *
import sys
import re
from itertools import izip
import numpy as np

parser = argparse.ArgumentParser(description='Analyze bulk lattices')
parser.add_argument('--models_re', '-m', action='store', type=str, help='models to include', default='*')
parser.add_argument('--tests_re', '-t', action='store', type=str, help='tests to include', default='*vacancy*,*interstitial*')
parser.add_argument('--label', '-l', action='store', help='optional label for models/tests directories', default='')
args = parser.parse_args()

try:
    with open("DEFAULTS_LABEL","r") as f:
        defaults_label = f.readline().strip()
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

multicomponent_constraints_data = get_multicomponent_constraints(args.label, models, default_analysis_settings["multicomponent_constraints"])

# print "multicomponent_constraints_data ", multicomponent_constraints_data

def mu_range(label, model_name, corresponding_bulk_test, multicomponent_constraints_data):
    (cur_min_EV, cur_composition) = read_ref_bulk_model_struct(label, model_name, corresponding_bulk_test)

    Leq = []
    Veq = cur_min_EV * len(cur_composition)
    print "equality constraint:",
    cur_elements = []
    i_of_element = {}
    for (i, element) in enumerate(cur_composition):
        if i > 0:
            print "+",
        print "{}*mu_{}".format(element[1], element[0]),
        cur_elements.append(element[0])
        Leq.append(element[1])
        i_of_element[element[0]] = i
    print " = mu_{} = {}".format(corresponding_bulk_test, cur_min_EV*len(cur_composition))

    ## print "constraint i_of_element, L, V ", i_of_element, Leq, "=", Veq

    Lne = []
    Vne = []
    print "inequalities:"
    for struct in multicomponent_constraints_data:
        # skip the current structure, since it is an equality constraint, not a stability limit inequality
        if struct == corresponding_bulk_test:
            continue
        # skip if any of the elements in this constraint aren't present in the material we're considering
        skip = False
        for element in multicomponent_constraints_data[struct]["composition"]:
            if not element[0] in cur_elements:
                skip = True
                continue
        if skip:
            continue

        # determine inequality
        Lne_cur = [0] * len(cur_composition)
        # print "   ", struct, multicomponent_constraints_data[struct]["E_per_atom"][model_name], multicomponent_constraints_data[struct]["composition"]
        n_atoms = 0
        for (i, element) in enumerate(multicomponent_constraints_data[struct]["composition"]):
            if i > 0:
                print "+",
            print "   {}*mu_{}".format(element[1], element[0]),
            n_atoms += element[1]
            Lne_cur[i_of_element[element[0]]] = element[1]
        print " <= mu_{} = {}".format(struct, multicomponent_constraints_data[struct]["E_per_atom"][model_name]*n_atoms)
        Lne.append(Lne_cur)
        Vne.append(multicomponent_constraints_data[struct]["E_per_atom"][model_name]*n_atoms)

    ## for (L, V) in izip(Lne, Vne):
        ## print "inequality i_of_element, L, V ", i_of_element, L, "<=", V

    mu_v_extrema = None
    if len(cur_composition) == 2: # binary
        Leq_3 = np.array([Leq[0], Leq[1], 0.0])
        mu_v_extrema = [ [ -np.finfo(float).max, None ] , [ np.finfo(float).max, None ]  ]
        for (L, V) in izip(Lne, Vne):
            A_lin_sys = np.array( [Leq, L] )
            b_lin_sys = np.array( [Veq, V] )
            x = np.linalg.solve(A_lin_sys, b_lin_sys)
            L_3 = np.array([L[0], L[1], 0.0])
            if np.cross(Leq_3, L_3)[2] < 0.0:
                # new max on component 0
                if x[0] < mu_v_extrema[1][0]:
                    mu_v_extrema[1][0] = x[0]
                    mu_v_extrema[1][1] = x[1]
            else:
                # new min on component 0
                if x[0] > mu_v_extrema[0][0]:
                    mu_v_extrema[0][0] = x[0]
                    mu_v_extrema[0][1] = x[1]

    elif len(cur_composition) == 3: # ternary
        raise("Can't do mu_range for %d-ary".format(len(cur_composition)))
    else:
        raise("Can't do mu_range for %d-ary".format(len(cur_composition)))

    if mu_v_extrema is None:
        mu_extrema = None
    else:
        mu_extrema = []
        for mu_pt in mu_v_extrema:
            mu_pt_Z = {}
            for Z in i_of_element:
                mu_pt_Z[Z] = mu_pt[i_of_element[Z]]
            mu_extrema.append(mu_pt_Z)

    return (mu_extrema)

# read and parse all data
data = {}
for model_name in models:
    sys.stderr.write("reading data for model {}\n".format(model_name))
    data[model_name] = {}
    cur_model_data = {}
    for test_name in tests:
        sys.stderr.write("   reading data for test {}\n".format(test_name))

        # read bulk test structure
        struct_filename_re = "{}model-{}-test-{}-ind_*_Z_*-relaxed.xyz".format(args.label, model_name, test_name)
        for struct_filename in glob.glob(struct_filename_re):
            try:
                struct = ase.io.read(struct_filename, format="extxyz")
            except:
                sys.stderr.write("No struct file '{}'\n".format(struct_filename))
                continue

            # read bulk test properties
            prop_filename ="{}model-{}-test-{}-properties.json".format(args.label, model_name, test_name)
            try:
                with open(prop_filename, "r") as model_data_file:
                    json_data = json.load(model_data_file)
            except:
                sys.stderr.write("No properties file '{}'\n".format(prop_filename))
                continue

            json_data["corresponding_bulk_test"] = struct.info["corresponding_bulk_test"]
            cur_model_data[test_name] = json_data.copy()

    data[model_name] = cur_model_data.copy()

for model_name in models:
    for test_name in tests:
        sys.stderr.write("do {} {}\n".format(model_name, test_name))

        mu_extrema = mu_range(args.label, model_name, data[model_name][test_name]["corresponding_bulk_test"], multicomponent_constraints_data)
        # print "mu_extrema", mu_extrema

        for vac in data[model_name][test_name]:
            if vac == "corresponding_bulk_test":
                continue
            print ""
            m = re.search("^ind_([0-9]*)_Z_([0-9]*)$", vac)
            ind = m.group(1)
            Z = m.group(2)
            if isinstance(data[model_name][test_name][vac], float) :
                print model_name, test_name, ind, Z, data[model_name][test_name][vac]
            else:
                Ef0 = data[model_name][test_name][vac][0]
                mu_str = data[model_name][test_name][vac][1]
                m = re.search('^\+mu_([0-9]*)$', mu_str)
                mu_Z = int(m.group(1))
                for mu_pt in mu_extrema:
                    print model_name, test_name, ind, Z, "multicomponent mu", mu_pt[mu_Z],"E_f", Ef0 + mu_pt[mu_Z]
