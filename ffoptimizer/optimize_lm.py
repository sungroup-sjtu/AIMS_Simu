#!/usr/bin/env python3
# coding=utf-8

import sys
import os
from collections import OrderedDict

import numpy as np
from lmfit import Parameters, Minimizer

sys.path.append('..')

from ffoptimizer.ppf import PPF
from ffoptimizer.db import DB
from ffoptimizer.target import Target, simulation

from sqlalchemy import and_


def build(ppf_file):
    for target in DB.session.query(Target).all():
        print(target.dir)
        target.build(ppf=ppf_file)


def optimize(ppf_file):
    def residual(params: Parameters):
        paras = OrderedDict()
        for k, v in params.items():
            paras[k] = v.value

        f = []
        for target in DB.session.query(Target).all():
            pres, hvap = target.get_pres_hvap_from_paras(ppf_file, paras)
            f.append((pres - target.P) * target.wDensity)
            f.append((hvap - target.hvap) * target.wHvap)

        return f

    def dfunc(params: Parameters):
        paras = OrderedDict()
        for k, v in params.items():
            paras[k] = v.value

        df = []
        for target in DB.session.query(Target).all():
            dPres, dHvap = target.get_dPres_dHvap_from_paras(ppf_file, paras)
            df.append([i * target.wDensity for i in dPres])
            df.append([i * target.wHvap for i in dHvap])
        print('\nDIRECTIVE:')
        for l in df:
            print(list(map(lambda x: round(x, 1), l)))

        return df

    def print_result(params: Parameters, iter: int, res: [float]):
        txt = '\nITERATION:%i\n' % iter
        txt += 'PARAMETERS:\n'
        for k, v in params.items():
            txt += '\t%s %10.5f\n' % (k, v.value)
        txt += 'RESIDUE: %s\n' % list(map(lambda x: round(x, 1), res))
        txt += 'RSQ: %.2f\n\n' % np.sum(list(map(lambda x: x ** 2, res)))

        print(txt)
        with open('Opt.log', 'a') as log:
            log.write(txt)

    params = Parameters()
    for k, v in adj_lj_paras.items():
        if k.endswith('r0'):
            params.add(k, value=v, min=2, max=5)
        elif k.endswith('e0'):
            params.add(k, value=v, min=0.01, max=3)
        elif k.endswith('bi'):
            params.add(k, value=v, min=-1, max=1)

    minimize = Minimizer(residual, params, iter_cb=print_result)
    result = minimize.leastsq(Dfun=dfunc, ftol=0.05)
    print(result.lmdif_message, '\n')

    return result.params


def run_npt(ppf_file):
    for target in DB.session.query(Target).all():
        target.run_npt(ppf_file)


def plot(ppfs):
    ppfs = [ppf[:-4] for ppf in ppfs]
    props = dict()
    for target in DB.session.query(Target).all():
        if not target.id in props.keys():
            props[target.id] = {'smiles': target.smiles, 'T': [], 'd_exp': [], 'h_exp': [],
                                'd_sim': OrderedDict(), 'h_sim': OrderedDict()}
        props[target.id]['T'].append(target.T)
        props[target.id]['d_exp'].append(target.density)
        props[target.id]['h_exp'].append(target.hvap)

        for ppf in ppfs:
            if ppf not in props[target.id]['d_sim'].keys():
                props[target.id]['d_sim'][ppf] = []
                props[target.id]['h_sim'][ppf] = []

            os.chdir(target.dir_npt)
            os.chdir(ppf)
            os.chdir('%i-%i' % (target.T, target.P))
            density = simulation.gmx.get_property('npt.edr', 'Density', begin=250)
            density /= 1000
            ei = simulation.gmx.get_property('hvap.edr', 'Potential', begin=250)
            hvap = 8.314 * target.T / 1000 - ei / target.n_mol

            props[target.id]['d_sim'][ppf].append(density)
            props[target.id]['h_sim'][ppf].append(hvap)

    os.chdir(sys.path[0])
    import pylab
    pylab.rcParams.update({'font.size': 16})
    for tid, prop in props.items():
        pylab.figure()
        pylab.plot(prop['T'], prop['d_exp'], '--')
        for ppf, points in prop['d_sim'].items():
            marker = 'x' if ppf.endswith('init') else 'o'
            pylab.plot(prop['T'], points, marker, label=ppf)
        y_mean = np.mean(prop['d_exp'])
        pylab.ylim(y_mean - 0.2, y_mean + 0.2)
        pylab.legend()
        pylab.title('Density %i %s (g/mL)' % (tid, prop['smiles']))
        pylab.savefig('density-%02i.png' % tid)

        pylab.figure()
        pylab.plot(prop['T'], prop['h_exp'], '--')
        for ppf, points in prop['h_sim'].items():
            marker = 'x' if ppf.endswith('init') else 'o'
            pylab.plot(prop['T'], points, marker, label=ppf)
        y_mean = np.mean(prop['h_exp'])
        pylab.ylim(y_mean - 10, y_mean + 10)
        pylab.legend()
        pylab.title('HVap %i %s (kJ/mol)' % (tid, prop['smiles']))
        pylab.savefig('hvap-%02i.png' % tid)


def init_db(filename):
    with open(filename) as f:
        lines = f.read().splitlines()
    for line in lines:
        if line.startswith('#') or line.strip() == '':
            continue
        words = line.split()

        target = Target()
        target.name = words[0]
        target.smiles = words[1]
        target.T = int(words[2])
        target.P = int(words[3])

        if DB.session.query(Target).filter(
                and_(Target.smiles == target.smiles,
                     Target.T == target.T,
                     Target.P == target.P
                     )).count() > 0:
            continue

        target.density = float(words[4])
        target.wDensity = float(words[5])
        target.hvap = float(words[6])
        target.wHvap = float(words[7])
        target.calc_n_mol()
        DB.session.add(target)

    try:
        DB.session.commit()
    except:
        DB.session.rollback()


if __name__ == '__main__':
    CWD = os.getcwd()
    DB.conn()
    cmd = sys.argv[1]
    if cmd == 'init':
        init_db(sys.argv[2])

    if cmd == 'build':
        ppf_file = None
        if len(sys.argv) == 3:
            ppf_file = os.path.abspath(sys.argv[2])
        build(ppf_file)

    if cmd == 'optimize':
        ppf_file = os.path.abspath(sys.argv[2])
        ppf = PPF(ppf_file)
        global adj_lj_paras
        adj_lj_paras = ppf.adj_lj_paras
        params_out = optimize(ppf_file)

        cycle = DB.session.query(Target).first().cycle + 1
        DB.session.query(Target).update({'cycle': cycle})
        DB.session.commit()

        paras = OrderedDict()
        for k, v in params_out.items():
            paras[k] = v.value
        ppf.set_lj_para(paras)
        ppf_out = os.path.join(CWD, 'TEAM-opt-%i.ppf' % cycle)
        ppf.write(ppf_out)
        build(os.path.abspath(ppf_out))

    if cmd == 'npt':
        ppf_file = os.path.abspath(sys.argv[2])
        run_npt(ppf_file)

    if cmd == 'plot':
        ppfs = sys.argv[2:]
        plot(ppfs)
