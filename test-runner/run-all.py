#!/usr/bin/env python

import os
import sys
import glob
import re
import json
import argparse

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

if args.label is not None:
    label_str = args.label
else:
    label_str=''

models = None
tests = None
omittests = None
if args.models is not None:
    models = [ os.path.join("../models/"+label_str, d) for d in args.models ]
if args.tests is not None:
    tests = [ os.path.join("../tests/"+label_str, d) for d in args.tests ]
if args.omit_tests is not None:
    omittests = [ os.path.join("../tests/"+label_str, d) for d in args.omit_tests ]
force = args.force
bugs = args.bugs

run_model_test = 'run-model-test.py'

if models is None:
    models = glob.glob(os.path.join('..', 'models/%s/*' % label_str))
if tests is None:
    tests = glob.glob(os.path.join('..', 'tests/%s/*' % label_str))

if omittests is not None:
    for d in omittests:
        try:
            tests.remove(d)
        except:
            sys.stderr.write("WARNING: omitted test %s not in list of tests\n" % d)
    
print 'Models:', models
print 'Tests:', tests

try:
    print "trying to read default_model_run_opts.json"
    with open("default_model_run_opts.json","r") as model_run_defaults_file:
        print "trying to parse json default_model_run_opts.json"
        model_run_defaults = json.load(model_run_defaults_file)
        print "got json ",model_run_defaults
except:
    model_run_defaults= {}

for model in models:

    model_name = os.path.split(model)[-1]

    print "Using model defaults ", model_run_defaults[model_name]
    try:
        use_bugs = 'use_bugs' in model_run_defaults[model_name]
    except:
        use_bugs = args.bugs
    try:
        use_mpi = 'use_mpi' in model_run_defaults[model_name]
    except:
        use_mpi = args.MPI
    try:
        use_openmp = 'use_openmp' in model_run_defaults[model_name]
    except:
        use_openmp = args.OpenMP

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
        if args.label is not None:
            cmd_args += ' --label '+args.label
        if args.base_model is not None:
            cmd_args += ' --base_model '+args.base_model

        if use_bugs:
            ident_string= '{0}_{1}'.format(model_name, test_name)
            if use_mpi:
                bugs_script='test.bugs_script_mpi'
                mpi_cmd=''
            else:
                mpi_cmd='-no_mpirun'
                bugs_script='test.bugs_script'
            if use_openmp:
                bugs_script+='_openmp'
            bugs_script += '.'+os.environ['HOSTNAME']
            cmd=('env REDIRECT_IO="'+cmd_args+'" bugs -exec=python '+mpi_cmd+' -time=96h -np='+('%d' % np)+' -script='+bugs_script+
                 ' -name='+ident_string+' -fileroot='+ident_string+' -in_cwd -output=job.'+ident_string)
        else:
            cmd="python "+cmd_args
        print cmd
        os.system(cmd)
