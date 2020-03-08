#!/usr/bin/env python3
# coding=utf-8
import sys
sys.path.append('..')
from app import create_app
from app.models import *
from app.models_cv import *
from app.selection import *
from mstools.formula import *
import argparse

parser = argparse.ArgumentParser(description='This is a code compare simulation data for stereoisomers')
parser.add_argument('-p', '--procedure', type=str, help='procedure of the compute: npt(ppm), or nvt-slab')
opt = parser.parse_args()

procedure = opt.procedure
app = create_app(procedure)
app.app_context().push()

if procedure == 'npt':
    f1 = open('density_compare.txt', 'w')
    f2 = open('cp_compare.txt', 'w')
    f3 = open('hvap_compare.txt', 'w')
    f4 = open('einter_compare.txt', 'w')
    f5 = open('compress_compare.txt', 'w')
    f6 = open('expansion_compare.txt', 'w')
elif procedure == 'nvt-slab':
    f1 = open('st_compare.txt', 'w')
    f2 = open('tc_compare.txt', 'w')
    f3 = open('dc_compare.txt', 'w')
else:
    sys.exit()
tasks = Task.query.filter(Task.procedure == procedure)
for task in tasks:
    if not task_selection(task):
        continue
    smiles = json.loads(task.smiles_list)[0]
    if has_stereo_isomer(smiles):
        stereo_smiles_list = get_stereo_isomer(smiles, canonical=True)
        py_mol = task.get_mol_list()[0]
        f = Formula(py_mol.formula)
        n_CON = f.atomdict.get('C', 0) + f.atomdict.get('O', 0) + f.atomdict.get('N', 0)
        cv = Cv.query.filter(Cv.smiles == smiles).first()
        t_list = json.loads(task.t_list)
        p_list = json.loads(task.p_list)
        for stereo_smiles in stereo_smiles_list:
            _task = tasks.filter(Task.smiles_list == json.dumps([stereo_smiles])).first()
            _cv = Cv.query.filter(Cv.smiles == stereo_smiles).first()
            if _task is not None:
                for i, t in enumerate(t_list):
                    if procedure == 'npt':
                        for p in p_list:
                            post_data = task.get_post_data(t, p)
                            density = post_data.get('density')
                            einter = post_data.get('einter')
                            hvap = (1 - n_CON / 15) * 8.314 * t / 1000 - einter if einter is not None else None
                            expansion = post_data.get('expansion')
                            compress = post_data.get('compress')
                            cp_inter = post_data.get('cp_inter')
                            cp_pv = post_data.get('cp_pv')
                            _post_data = _task.get_post_data(t, p)
                            _density = _post_data.get('density')
                            _einter = _post_data.get('einter')
                            _hvap = (1 - n_CON / 15) * 8.314 * t / 1000 - _einter if _einter is not None else None
                            _expansion = _post_data.get('expansion')
                            _compress = _post_data.get('compress')
                            _cp_inter = _post_data.get('cp_inter')
                            _cp_pv = _post_data.get('cp_pv')
                            if None not in [density, _density]:
                                f1.write('%f %f %s %f % s %f\n' % (t, p, smiles, density, stereo_smiles, _density))
                            if None not in [cv, cp_inter, cp_pv, _cv, _cp_inter, _cp_pv]:
                                cp = cv.get_post_cv(t) + cp_inter + cp_pv
                                _cp = _cv.get_post_cv(t) + _cp_inter + _cp_pv
                                f2.write('%f %f %s %f % s %f\n' % (t, p, smiles, cp, stereo_smiles, _cp))
                            if None not in [hvap, _hvap]:
                                f3.write('%f %f %s %f % s %f\n' % (t, p, smiles, hvap, stereo_smiles, _hvap))
                            if None not in [einter, _einter]:
                                f4.write('%f %f %s %f % s %f\n' % (t, p, smiles, einter, stereo_smiles, _einter))
                            if None not in [compress, _compress]:
                                f5.write('%f %f %s %f % s %f\n' % (t, p, smiles, compress, stereo_smiles, _compress))
                            if None not in [expansion, _expansion]:
                                f6.write('%f %f %s %f % s %f\n' % (t, p, smiles, expansion, stereo_smiles, _expansion))
                    elif procedure == 'nvt-slab':
                        post_data = task.get_post_data(t)
                        st = post_data.get('st')
                        _post_data = _task.get_post_data(t)
                        _st = _post_data.get('st')
                        if None not in [st, _st]:
                            f1.write('%f %s %f % s %f\n' % (t, smiles, st, stereo_smiles, _st))
                        if i == 0:
                            tc = post_data.get('tc')
                            dc = post_data.get('dc')
                            _tc = _post_data.get('tc')
                            _dc = _post_data.get('dc')
                            if None not in [tc, _tc]:
                                f2.write('%s %f % s %f\n' % (smiles, tc, stereo_smiles, _tc))
                            if None not in [dc, _dc]:
                                f3.write('%s %f % s %f\n' % (smiles, dc, stereo_smiles, _dc))