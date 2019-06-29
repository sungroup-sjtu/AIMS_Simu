#!/usr/bin/env python3
# coding=utf-8

import sys
sys.path.append('..')
from app import create_app

import argparse
parser = argparse.ArgumentParser(description='This is a code to post-process the results related to specific compute')
parser.add_argument('-p', '--procedure', type=str, help='procedure of the compute: npt(ppm), or nvt-slab')
parser.add_argument('-o', '--overwrite', type=bool, help='overwrite the existent post-result', default=False)
opt = parser.parse_args()

procedure = opt.procedure
app = create_app(procedure)
app.app_context().push()

from app.models import *
tasks = Task.query.filter(Task.procedure == procedure)
print('Total task = %i' % (tasks.count()))
count = 1
for task in tasks:
    sys.stdout.write('\rpresent task = %i\n' % (count))
    count += 1
    task.post_process(overwrite=opt.overwrite)
    task.get_LJ_atom_type()
