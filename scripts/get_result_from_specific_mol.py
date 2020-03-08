#!/usr/bin/env python3
# coding=utf-8

import sys
sys.path.append('..')
from app import create_app
import argparse
parser = argparse.ArgumentParser(description='This is a code to post-process the results related to specific compute')
parser.add_argument('-p', '--procedure', type=str, help='procedure of the compute: npt(ppm), or nvt-slab')
parser.add_argument('-s', '--smiles', type=str, help='The smiles list of specific molecule')
parser.add_argument('-pro', '--property', type=str, help='The property want to calculated')
opt = parser.parse_args()

app = create_app(opt.procedure)
app.app_context().push()

from app.models import *
smiles_list = opt.smiles.split('.')
# task_id = Task.query.filter(Task.smiles_list == json.dumps(smiles_list)).first().id
print('# T P sim simstderr')
for job in Job.query.join(Task).filter(Task.smiles_list == json.dumps(smiles_list)).order_by(Job.t, Job.p):
# for job in Job.query.filter(Job.task_id == task_id).order_by(Job.t, Job.p):
    if job.converged:
        if json.loads(job.result).get(opt.property) != None:
            value, v_stderr = json.loads(job.result).get(opt.property)
            print('%.2f %.1f %f %f' % (job.t, job.p, value, v_stderr))
