#!/usr/bin/env python

from __future__ import print_function


import os
import sys
import datetime
import json
import time
import argparse

if 'USE_MPI' in os.environ:
    import mpi4py
    from quippy.mpi_context import MPI_context
    do_io = (mpi4py.MPI.COMM_WORLD.Get_rank() == 0)
    print("rank ",mpi4py.MPI.COMM_WORLD.Get_rank(), "io", do_io)
else:
    do_io = True

try:
    from ase.calculators.checkpoint import CheckpointCalculator
except:
    pass

#import logging
#logging.basicConfig(stream=sys.stdout, level=logging.INFO)
#logger = logging.getLogger('ase.optimize.precon')
#logger.propagate = True
#logger.setLevel(logging.DEBUG)


parser = argparse.ArgumentParser(description='Run a particular model-test combination')
parser.add_argument('model', type=str, action='store', help='model name')
parser.add_argument('test', type=str, action='store', help='test name')
parser.add_argument('--system_label','-l', type=str, action='store', help='label for tests directory', default='')
parser.add_argument('--force','-f', action='store_true', help='force rerunning of test')
parser.add_argument('--base_model','-B', type=str, action='store', help='optional base model to start from')
parser.add_argument('--no_redirect_io','-N', action='store_true', help='do not redirect io')
parser.add_argument('--no_append_log', dest='append_log', action='store_false', help='overwrite log (.txt) file, rather than appending')
args = parser.parse_args()

model_path = os.path.split(args.model)[0]
if len(model_path) == 0:
    model_path = os.path.join("..","..","models")
model_name = os.path.split(args.model)[1]

test_path = os.path.split(args.test)[0]
if len(test_path) == 0:
    test_path = os.path.join("..","..","tests")
test_name = os.path.split(args.test)[1]

force = args.force

my_path = os.path.split(os.path.realpath(__file__))[0]
# read utilities from current relative path
share_dir = os.path.join(my_path, '..', 'share')
print("share_dir",share_dir)
sys.path.insert(0, share_dir)
import utilities
## # remove this path, since later relative path will be different
## sys.path.remove(share_dir)

# set quantities to be available to utilities routines
utilities.model_name = model_name
utilities.base_model_name = args.base_model
utilities.test_name = test_name
utilities.system_label = args.system_label

run_root = utilities.model_test_root()
dir_name = "run_"+run_root # maybe a different name?
if not os.path.exists(dir_name) and do_io:
    os.mkdir(dir_name)
time.sleep(1)
os.chdir(dir_name)

utilities.run_root = run_root
utilities.base_run_root = utilities.model_test_root(base_model=True)

json_file_name = os.path.join('..', run_root+'-properties.json')
if not force and os.path.isfile(json_file_name) and os.path.getsize(json_file_name) > 0:
    print("%s already exists and is not empty, not rerunning test" % json_file_name)
    sys.exit(0)

model_dir = os.path.join(model_path, model_name)
test_dir = os.path.join(test_path, test_name)

## sys.path.insert(0, share_dir)
sys.path.insert(0, model_dir)
sys.path.insert(0, test_dir)

if not args.no_redirect_io:
    _stdout, _stderr = sys.stdout, sys.stderr
    if do_io:
        if args.append_log:
            log = open(os.path.join('..',run_root+'.txt'), 'a', 1)
        else:
            log = open(os.path.join('..',run_root+'.txt'), 'w', 1)
        sys.stdout, sys.stderr = log, log
    else:
        sys.stdout = open(os.devnull, "w")
        sys.stderr = open(os.devnull, "w")
sys.stdout.write('\n###### START RUN '+time.ctime()+' ######\n')

import logging
logging.basicConfig(stream=sys.stdout, level=logging.INFO)
logger = logging.getLogger('ase.optimize.precon')
logger.propagate = True
logger.setLevel(logging.INFO)

print('Model {0}, Test {1}'.format(model_name, test_name))
print('Test run at {0}'.format(datetime.datetime.now().strftime("%Y-%m-%d %H:%M")))
print('')

model_file = os.path.join(model_dir, 'model.py')
print('model file:',model_file)
print('='*60)
sys.stdout.write(open(model_file).read())
print('='*60)

test_file = os.path.join(test_dir, 'test.py')
print('test file:', test_file)
print('='*60)
sys.stdout.write(open(test_file).read())
print('='*60)

import model # import and run the current model

# adapt model to test, e.g. to set output directory name
if hasattr(model, 'start'):
    model.start(test_name)

# create a checkpoint file for this specific (model, test) combination
if ('GAP_TESTER_CHECKPOINT' in os.environ and os.environ['GAP_TESTER_CHECKPOINT'] != "") or not hasattr(model, 'no_checkpoint') or not model.no_checkpoint:
   checkpoint_file = run_root+'.db'
   print('Using checkpoint file',checkpoint_file)
   model.calculator = CheckpointCalculator(model.calculator, db=checkpoint_file)

import test  # import and run the current test

print('='*60)
print('Property calculation output:')
print('')

# serialise results in machine readable JSON format
json_file = open(json_file_name, 'w')
json.dump(test.properties, json_file)
json_file.close()

print('')
print('Summary of computed properties:')
print(test.properties)

print('='*60)

if not args.no_redirect_io:
    sys.stdout, sys.stderr = _stdout, _stderr
    if do_io:
        log.close()

if hasattr(model, 'shutdown'):
    model.shutdown()
