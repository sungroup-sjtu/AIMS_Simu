#!/usr/bin/env python3
# coding=utf-8

import sys

sys.path.append('..')
from app import create_app
from app.models import *
from app.models_nist import *
from app.models_cv import Cv
from app.selection import *
import argparse

parser = argparse.ArgumentParser(description='This is a code to generate input txt from experimental data')
parser.add_argument('-t', '--type', type=str, help='The type of the output(SIM, EXP)')
parser.add_argument('--database', type=str, help='The type of the experimental database(ILTHERMO, NIST)')
parser.add_argument('--procedure', type=str, help='The procedure of simulation(npt)')
parser.add_argument('--selection', type=bool, help='Use the selection function', default=True)
parser.add_argument('--nheavy', type=int, help='Only 5 < heavy atoms < nheavy is considered', default=16)
parser.add_argument('--training', type=str, help='The txt file for training sets smiles', default=None)
parser.add_argument('--errormolecules', type=bool, help='Find the molecules in training set without simulation data',
                    default=False)
parser.add_argument('--rawdata', type=bool, help='Use the raw simulation data, otherwise fit data', default=False)
parser.add_argument('--explicit', type=bool, help='Only use the simulation data with unambiguous smiles', default=True)

opt = parser.parse_args()
app = create_app(opt.procedure)
app.app_context().push()


def write_data(_smiles, _training_smiles_list, _t, _p, _property, file, file_train, positive=False):
    if _property is not None:
        if positive and _property < 0:
            return
        if _t is None:
            file.write('%s %.3e\n' % (smiles, _property))
            if _smiles in _training_smiles_list:
                file_train.write('%s %.3e\n' % (smiles, _property))
        elif _p is None:
            file.write('%s %i %.3e\n' % (smiles, _t, _property))
            if _smiles in _training_smiles_list:
                file_train.write('%s %i %.3e\n' % (smiles, _t, _property))
        else:
            file.write('%s %i %i %.3e\n' % (smiles, _t, _p, _property))
            if _smiles in _training_smiles_list:
                file_train.write('%s %i %i %.3e\n' % (smiles, _t, _p, _property))


if opt.type == 'EXP' and opt.database == 'NIST':
    f1 = open('Tvap.txt', 'w')
    f1.write('SMILES Tvap\n')
    f2 = open('Tm.txt', 'w')
    f2.write('SMILES Tm\n')
    f3 = open('Tc.txt', 'w')
    f3.write('SMILES Tc\n')
    molecules = NistMolecule.query.filter(NistMolecule.n_heavy > 5)
    for mol in molecules:
        print(mol.smiles)
        if not selection(mol.smiles, type='CH'):
            continue
        if mol.tb is not None:
            f1.write('%s %f\n' % (mol.smiles, mol.tb))
        if mol.tt is not None:
            f2.write('%s %f\n' % (mol.smiles, mol.tt))
        if mol.tc is not None:
            f3.write('%s %f\n' % (mol.smiles, mol.tc))
elif opt.type == 'SIM' and opt.errormolecules and opt.training:
    info = pd.read_csv(opt.training, sep='\s+', header=0)
    training_smiles_list = []
    for smiles in info['SMILES'].unique():
        training_smiles_list.append(get_canonical_smiles(smiles))
    print('%i molecules in training set' % len(training_smiles_list))
    for smiles in training_smiles_list:
        task = Task.query.filter(Task.procedure == opt.procedure).filter(
            Task.smiles_list == json.dumps([smiles])).first()
        if task is None:
            print('There is no task for %s' % smiles)
            continue
        if task.status in [Compute.Status.FAILED, Compute.Status.STARTED]:
            print('Task for %s not finished' % smiles)
            continue
        if not task_selection(task, select=opt.selection):
            print('%s does not satisfy selection rule' % smiles)
            continue
        post_result = task.get_post_result()
        if post_result is None:
            print('There is no post_result for %s' % smiles)
            continue
        t = task.get_t_list()[0]
        post_data = task.get_post_data(T=t)
        if opt.procedure == 'npt':
            density = post_data.get('density')
            einter = post_data.get('einter')
            expansion = post_data.get('expansion')
            cp_inter = post_data.get('cp_inter')
            cp_pv = post_data.get('cp_pv')
            cv = Cv.query.filter(Cv.smiles == smiles).first()
            if density is None:
                print('%s has no density' % json.loads(task.smiles_list))
            if einter is None:
                print('%s has no hvap' % json.loads(task.smiles_list))
            if expansion is None:
                print('%s has no expansion' % json.loads(task.smiles_list))
            if cv is None:
                print('%s has no dft cv' % json.loads(task.smiles_list))
            if None in [cp_inter, cp_pv]:
                print('%s has no md cv' % json.loads(task.smiles_list))
        elif opt.procedure == 'nvt-slab':
            st = post_data.get('st')
            tc = post_data.get('tc')
            dc = post_data.get('dc')
            if st is None:
                print('%s has no st' % json.loads(task.smiles_list))
            if tc is None:
                print('%s has no tc' % json.loads(task.smiles_list))
            if dc is None:
                print('%s has no dc' % json.loads(task.smiles_list))
elif opt.type == 'SIM' and opt.procedure == 'npt':
    if opt.training is not None:
        f_train1 = open('result-ML-density-train.txt', 'w')
        f_train2 = open('result-ML-einter-train.txt', 'w')
        f_train3 = open('result-ML-compress-train.txt', 'w')
        f_train4 = open('result-ML-hvap-train.txt', 'w')
        f_train5 = open('result-ML-cp-train.txt', 'w')
        f_train6 = open('result-ML-expansion-train.txt', 'w')
        f_train1.write('SMILES T P density\n')
        f_train2.write('SMILES T P einter\n')
        f_train3.write('SMILES T P compress\n')
        f_train4.write('SMILES T P hvap\n')
        f_train5.write('SMILES T P cp\n')
        f_train6.write('SMILES T P expansion\n')
        info = pd.read_csv(opt.training, sep='\s+', header=0)
        training_smiles_list = []
        for smiles in info['SMILES'].to_list():
            training_smiles_list.append(get_canonical_smiles(smiles))
    else:
        training_smiles_list = []
        f_train1 = f_train2 = f_train3 = f_train4 = f_train5 = f_train6 = None
    f1 = open('result-ML-density.txt', 'w')
    f2 = open('result-ML-einter.txt', 'w')
    f3 = open('result-ML-compress.txt', 'w')
    f4 = open('result-ML-hvap.txt', 'w')
    f5 = open('result-ML-cp.txt', 'w')
    f6 = open('result-ML-expansion.txt', 'w')
    f1.write('SMILES T P density\n')
    f2.write('SMILES T P einter\n')
    f3.write('SMILES T P compress\n')
    f4.write('SMILES T P hvap\n')
    f5.write('SMILES T P cp\n')
    f6.write('SMILES T P expansion\n')
    smiles_list = []
    t_list = []
    p_list = []
    den_list = []
    tasks = Task.query.filter(Task.procedure == opt.procedure)
    for task in tasks:
        if task.status in [Compute.Status.FAILED, Compute.Status.STARTED]:
            continue
        if not task_selection(task, select=opt.selection):
            continue
        if opt.explicit and has_stereo_isomer(json.loads(task.smiles_list)[0]):
            continue
        post_result = task.get_post_result()
        if post_result is None or post_result['density-poly4'][-1] < 0.999:
            continue
        smiles = task.get_smiles_list()[0]
        py_mol = task.get_mol_list()[0]
        if not 5 < get_heavy_atom_numbers(smiles) < opt.nheavy:
            continue
        f = Formula(py_mol.formula)
        n_CON = f.atomdict.get('C', 0) + f.atomdict.get('O', 0) + f.atomdict.get('N', 0)
        cv = Cv.query.filter(Cv.smiles == smiles).first()
        jobs = task.jobs
        n_mol_list = json.loads(task.n_mol_list)
        for t in task.get_t_list():
            for p in task.get_p_list():
                job_result = json.loads(jobs.filter(Job.t == t).filter(Job.p == p).first().result)
                post_data = task.get_post_data(t, p)
                if opt.rawdata and job_result.get('density') is not None:
                    density = job_result.get('density')[0]
                else:
                    density = post_data.get('density')
                if opt.rawdata and job_result.get('einter') is not None:
                    einter = job_result.get('einter')[0] / min(n_mol_list)
                else:
                    einter = post_data.get('einter')
                if opt.rawdata and job_result.get('expansion') is not None:
                    expansion = job_result.get('expansion')[0]
                else:
                    expansion = post_data.get('expansion')
                if opt.rawdata and job_result.get('compress') is not None:
                    compress = job_result.get('compress')[0]
                else:
                    compress = post_data.get('compress')

                if einter is not None:
                    hvap = (1 - n_CON / 15) * 8.314 * t / 1000 - einter
                else:
                    hvap = None

                write_data(smiles, training_smiles_list, t, p, density, f1, f_train1, positive=True)
                write_data(smiles, training_smiles_list, t, p, einter, f2, f_train2)
                write_data(smiles, training_smiles_list, t, p, compress, f3, f_train3, positive=True)
                write_data(smiles, training_smiles_list, t, p, hvap, f4, f_train4, positive=True)
                write_data(smiles, training_smiles_list, t, p, expansion, f6, f_train6, positive=True)

                cp_inter = post_data.get('cp_inter')
                cp_pv = post_data.get('cp_pv')
                if None not in [cv, cp_inter, cp_pv]:
                    cv_intra = cv.get_post_cv(t)
                    cp = cv_intra + cp_inter + cp_pv
                    write_data(smiles, training_smiles_list, t, p, cp, f5, f_train5, positive=True)
elif opt.type == 'SIM' and opt.procedure == 'nvt-slab':
    if opt.training is not None:
        f_train1 = open('result-ML-st-train.txt', 'w')
        f_train2 = open('result-ML-tc-train.txt', 'w')
        f_train3 = open('result-ML-dc-train.txt', 'w')
        f_train1.write('SMILES T st\n')
        f_train2.write('SMILES tc\n')
        f_train3.write('SMILES dc\n')
        info = pd.read_csv(opt.training, sep='\s+', header=0)
        training_smiles_list = []
        for smiles in info['SMILES'].to_list():
            training_smiles_list.append(get_canonical_smiles(smiles))
    else:
        f_train1 = f_train2 = f_train3 = None
        training_smiles_list = []
    f1 = open('result-ML-st.txt', 'w')
    f2 = open('result-ML-tc.txt', 'w')
    f3 = open('result-ML-dc.txt', 'w')
    f1.write('SMILES T st\n')
    f2.write('SMILES tc\n')
    f3.write('SMILES dc\n')
    t_list = []
    p_list = []
    tasks = Task.query.filter(Task.procedure == opt.procedure)
    for task in tasks:
        if task.status in [Compute.Status.FAILED, Compute.Status.STARTED]:
            continue
        if not task_selection(task, select=opt.selection):
            continue
        if opt.explicit and has_stereo_isomer(json.loads(task.smiles_list)[0]):
            continue
        post_result = task.get_post_result()
        if post_result is None:
            continue
        smiles = task.get_smiles_list()[0]
        if not 5 < get_heavy_atom_numbers(smiles) < opt.nheavy:
            continue
        jobs = task.jobs
        for i, t in enumerate(task.get_t_list()):
            job_result = json.loads(jobs.filter(Job.t == t).first().result)
            post_data = task.get_post_data(T=t)
            if opt.rawdata and job_result.get('st') is not None:
                st = job_result.get('st')[0]
            else:
                st = post_data.get('st')

            write_data(smiles, training_smiles_list, t, None, st, f1, f_train1, positive=True)

            if i == 0:
                tc = post_data.get('tc')
                write_data(smiles, training_smiles_list, None, None, tc, f2, f_train2, positive=True)
                dc = post_data.get('dc')
                write_data(smiles, training_smiles_list, None, None, dc, f3, f_train3, positive=True)
