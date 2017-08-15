#!/usr/bin/env python

import os
import sys
import glob
import re
import json
import argparse
import itertools

parser = argparse.ArgumentParser(description='Run all tests on all models')
parser.add_argument('--models', '-m', action='store', nargs='*', type=str, help='models to include')
parser.add_argument('--tests', '-t', action='store', nargs='*', type=str, help='tests to include')
parser.add_argument('--omit-tests', '-o', action='store', nargs='*', type=str, help='tests to omit')
parser.add_argument('--force', '-f', action='store_true', help='force rerunning of tests')
parser.add_argument('--bugs', '-b', action='store_true', help='use bugs to generate job scripts')
parser.add_argument('--MPI', '-M', action='store_true', help='use MPI')
parser.add_argument('--OpenMP', '-O', action='store_true', help='use OpenMP')
parser.add_argument('--label', '-l', action='store', help='optional label for models/tests directories')
parser.add_argument('--base_model', '-B', action='store', type=str, help='model to use as initial config for tests where it is enabled')
args = parser.parse_args()


try:
    with open("DEFAULTS_LABEL","r") as f:
        defaults_label = f.readline().strip()+"_"
except:
    defaults_label = ""

try:
    print "reading","{}default_run_opts.json".format(defaults_label)
    with open("{}default_run_opts.json".format(defaults_label)) as defaults_file:
        defaults_vals = json.load(defaults_file)
        if 'models' in defaults_vals and args.models is None:
            args.models = defaults_vals['models']
        if 'tests' in defaults_vals and args.tests is None:
            args.tests = defaults_vals['tests']
        if 'omit-tests' in defaults_vals and args.omit_tests is None:
            args.omit_tests = defaults_vals['omit_tests']
        if 'label' in defaults_vals and args.label is None:
            args.label = defaults_vals['label']
        if 'base_model' in defaults_vals and args.base_model is None:
            args.base_model = defaults_vals['base_model']
        if 'bugs' in defaults_vals:
            args.bugs = defaults_vals['bugs']
except:
    sys.stderr.write("No parsable {}default_run_opts.json file\n".format(defaults_label))

if args.label is None:
    args.label = ''

models = None
tests = None
omittests = None
if args.models is not None:
    models = list(itertools.chain.from_iterable([ glob.glob(os.path.join("../models/"+args.label, d)) for d in args.models ]))
if args.tests is not None:
    tests = list(itertools.chain.from_iterable([ glob.glob(os.path.join("../tests/"+args.label, d)) for d in args.tests ]))
if args.omit_tests is not None:
    omittests = list(itertools.chain.from_iterable([ glob.glob(os.path.join("../tests/"+args.label, d)) for d in args.omit_tests ]))
force = args.force
bugs = args.bugs

run_model_test = 'run-model-test.py'

if models is None:
    models = glob.glob(os.path.join('..', 'models/%s/*' % args.label))
if tests is None:
    tests = glob.glob(os.path.join('..', 'tests/%s/*' % args.label))

if omittests is not None:
    for d in omittests:
        try:
            tests.remove(d)
        except:
            sys.stderr.write("WARNING: omitted test %s not in list of tests\n" % d)
    
print 'Models:', models
print 'Tests:', tests

try:
    print "trying to read {}default_model_run_opts.json".format(defaults_label)
    with open("{}default_model_run_opts.json".format(defaults_label),"r") as model_run_defaults_file:
        model_run_defaults = json.load(model_run_defaults_file)
except:
    model_run_defaults= {}
print "got model_run_defaults", model_run_defaults

for model in models:

    model_name = os.path.split(model)[-1]

    use_model_run_defaults = {}
    # check for exact matches
    for model_default_key in model_run_defaults:
        if model_name == model_default_key:
            use_model_run_defaults = model_run_defaults[model_default_key]
    # check for regexp
    if len(use_model_run_defaults) == 0:
        for model_default_key in model_run_defaults:
            if re.search(model_default_key, model_name) is not None:
                use_model_run_defaults = model_run_defaults[model_default_key]
                break

    if not args.MPI:
        try:
            args.MPI = 'MPI' in use_model_run_defaults
        except:
            pass

    if not args.OpenMP:
        try:
            args.OpenMP = 'OpenMP' in use_model_run_defaults
        except:
            pass

    for test in tests:
        test_name = os.path.split(test)[-1]
        try:
            fp = open(model+"/COMPUTATIONAL_COST")
            model_cost = float(fp.readline().strip())
        except:
            model_cost = 1.0
        try:
            fp = open(test+"/COMPUTATIONAL_COST")
            test_cost = float(fp.readline().strip())
        except:
            test_cost = 1.0
        np=int(round(int(test_cost*model_cost*16)/16.0))*16
        if np < 16:
            np = 16

        cmd_args = '{0} {1} {2}'.format(run_model_test, model_name, test_name)
        if force:
            cmd_args += ' -force'
        if args.label is not None and len(args.label) > 0:
            cmd_args += ' --label '+args.label
        if args.base_model is not None:
            cmd_args += ' --base_model '+args.base_model

        if args.bugs:
            if args.label is not None:
                ident_string= '{0}_{1}_{2}'.format(args.label, model_name, test_name)
            else:
                ident_string= '{0}_{1}'.format(model_name, test_name)
            if args.MPI:
                bugs_script='test.bugs_script_mpi'
                mpi_cmd=''
            else:
                mpi_cmd='-no_mpirun'
                bugs_script='test.bugs_script'
            if args.OpenMP:
                bugs_script+='_openmp'
            bugs_script += '.'+os.environ['HOSTNAME']
            cmd=('env REDIRECT_IO="'+cmd_args+'" bugs -exec=python '+mpi_cmd+' -time=96h -np='+('%d' % np)+' -script='+bugs_script+
                 ' -name='+ident_string+' -fileroot='+ident_string+' -in_cwd -output=job.'+ident_string)
        else:
            cmd="python "+cmd_args
        print cmd
        os.system(cmd)
