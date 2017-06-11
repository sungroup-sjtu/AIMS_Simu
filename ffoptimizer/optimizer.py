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
from ffoptimizer.target import Target

from sqlalchemy import and_

CWD = os.getcwd()


class Optimizer():
    def __init__(self, db_file):
        self.db = DB(db_file)
        self.db.conn()

    def build_nvt(self, ppf_file):
        for target in self.db.session.query(Target).all():
            print(target.dir_nvt)
            target.build_nvt(ppf=ppf_file)

    def optimize_nvt(self, ppf_file):
        def residual(params: Parameters):
            paras = OrderedDict()
            for k, v in params.items():
                paras[k] = v.value

            f = []
            for target in self.db.session.query(Target).all():
                pres, hvap = target.get_pres_hvap_from_paras(ppf_file, paras)
                f.append((pres - target.P) * target.wDensity)
                f.append((hvap - target.hvap) * target.wHvap)

            return f

        def dfunc(params: Parameters):
            paras = OrderedDict()
            for k, v in params.items():
                paras[k] = v.value

            df = []
            for target in self.db.session.query(Target).all():
                dPres, dHvap = target.get_dPres_dHvap_from_paras(ppf_file, paras)
                df.append([i * target.wDensity for i in dPres])
                df.append([i * target.wHvap for i in dHvap])

            txt = '\nDERIVATIVE:\n'
            for l in df:
                txt += '%s\n' % list(map(lambda x: round(x, 1), l))

            print(txt)
            with open(os.path.join(CWD, 'Opt.log'), 'a') as log:
                log.write(txt)

            return df

        def print_result(params: Parameters, iter: int, res: [float]):
            txt = '\nITERATION:%i\n' % iter
            txt += 'PARAMETERS:\n'
            for k, v in params.items():
                txt += '\t%s %10.5f\n' % (k, v.value)
            txt += 'RESIDUE: %s\n' % list(map(lambda x: round(x, 1), res))
            txt += 'RSQ: %.2f\n\n' % np.sum(list(map(lambda x: x ** 2, res)))

            print(txt)
            with open(os.path.join(CWD, 'Opt.log'), 'a') as log:
                log.write(txt)

        ppf = PPF(ppf_file)
        params = Parameters()
        for k, v in ppf.adj_lj_paras.items():
            if k.endswith('r0'):
                params.add(k, value=v, min=2, max=5)
            elif k.endswith('e0'):
                params.add(k, value=v, min=0.01, max=3)
            elif k.endswith('bi'):
                params.add(k, value=v, min=-1, max=1)

        minimize = Minimizer(residual, params, iter_cb=print_result)
        result = minimize.leastsq(Dfun=dfunc, ftol=0.1)
        print(result.lmdif_message, '\n')

        return result.params

    def optimize_npt(self, ppf_file):
        def residual(params: Parameters):
            paras = OrderedDict()
            for k, v in params.items():
                paras[k] = v.value

            res = []
            for target in self.db.session.query(Target).all():
                dens, hvap = target.get_npt_result(os.path.basename(ppf_file)[:-4])
                res.append((dens - target.density) / target.density * 100 * target.wDensity)
                res.append((hvap - target.hvap) / target.hvap * 100 * target.wHvap)

            return res

        def derivative(params: Parameters):
            paras = OrderedDict()
            for k, v in params.items():
                paras[k] = v.value

            jacobian = []
            for target in self.db.session.query(Target).all():
                dDens, dHvap = target.get_dDens_dHvap_from_paras(ppf_file, paras)
                jacobian.append([i / target.density * 100 * target.wDensity for i in dDens])
                jacobian.append([i / target.hvap * 100 * target.wHvap for i in dHvap])

            return jacobian

        ppf = PPF(ppf_file)
        params = Parameters()
        for k, v in ppf.adj_lj_paras.items():
            if k.endswith('r0'):
                params.add(k, value=v, min=2, max=5)
            elif k.endswith('e0'):
                params.add(k, value=v, min=0.01, max=3)
            elif k.endswith('bi'):
                params.add(k, value=v, min=-1, max=1)

        res = residual(params)
        jacobian = derivative(params)

        txt = '\nResidual:\n'
        for r in res:
            txt += '%f\n' % r
        txt = '\nJacobian Matrix:\n'
        for l in jacobian:
            txt += '%s\n' % l
        print(txt)
        with open(os.path.join(CWD, 'Opt.log'), 'a') as log:
            log.write(txt)

    def run_npt(self, ppf_file):
        for target in self.db.session.query(Target).all():
            target.run_npt(ppf_file)

    def plot(self, ppfs):
        ppfs = [ppf[:-4] for ppf in ppfs]
        props = dict()
        for target in self.db.session.query(Target).all():
            if not target.name in props.keys():
                props[target.name] = {'smiles': target.smiles, 'T': [], 'd_exp': [], 'h_exp': [],
                                      'd_sim': OrderedDict(), 'h_sim': OrderedDict()}
            props[target.name]['T'].append(target.T)
            props[target.name]['d_exp'].append(target.density)
            props[target.name]['h_exp'].append(target.hvap)

            for ppf in ppfs:
                if ppf not in props[target.name]['d_sim'].keys():
                    props[target.name]['d_sim'][ppf] = []
                    props[target.name]['h_sim'][ppf] = []

                density, hvap = target.get_npt_result(ppf)

                props[target.name]['d_sim'][ppf].append(density)
                props[target.name]['h_sim'][ppf].append(hvap)

        os.chdir(CWD)
        import pylab
        pylab.rcParams.update({'font.size': 12})
        for tid, prop in props.items():
            pylab.figure(figsize=(6, 8))
            pylab.subplot(211)
            pylab.plot(prop['T'], prop['d_exp'], '--')
            for ppf, points in prop['d_sim'].items():
                marker = 'x' if ppf.endswith('init') else 'o'
                pylab.plot(prop['T'], points, marker, label=ppf)
            y_mean = np.mean(prop['d_exp'])
            pylab.ylim(y_mean - 0.2, y_mean + 0.2)
            pylab.legend()
            pylab.title('Density %s %s (g/mL)' % (tid, prop['smiles']))

            pylab.subplot(212)
            pylab.plot(prop['T'], prop['h_exp'], '--')
            for ppf, points in prop['h_sim'].items():
                marker = 'x' if ppf.endswith('init') else 'o'
                pylab.plot(prop['T'], points, marker, label=ppf)
            y_mean = np.mean(prop['h_exp'])
            pylab.ylim(y_mean - 20, y_mean + 20)
            pylab.legend()
            pylab.title('HVap %s %s (kJ/mol)' % (tid, prop['smiles']))
            pylab.savefig('property-%s.png' % tid)

    def init_db(self, filename):
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

            if self.db.session.query(Target).filter(
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
            self.db.session.add(target)

        try:
            self.db.session.commit()
        except:
            self.db.session.rollback()


if __name__ == '__main__':
    optimizer = Optimizer(db_file='ffoptimizer.db')

    cmd = sys.argv[1]
    if cmd == 'init':
        optimizer.init_db(sys.argv[2])

    if cmd == 'build-nvt':
        ppf_file = None
        if len(sys.argv) == 3:
            ppf_file = os.path.abspath(sys.argv[2])
        optimizer.build_nvt(ppf_file)

    if cmd == 'optimize-nvt':
        ppf_file = os.path.abspath(sys.argv[2])
        params_out = optimizer.optimize_nvt(ppf_file)

        cycle = optimizer.db.session.query(Target).first().cycle + 1
        optimizer.db.session.query(Target).update({'cycle': cycle})
        optimizer.db.session.commit()

        paras = OrderedDict()
        for k, v in params_out.items():
            paras[k] = v.value

        ppf = PPF(ppf_file)
        ppf.set_lj_para(paras)
        ppf_out = os.path.join(CWD, 'TEAM-opt-%i.ppf' % cycle)
        ppf.write(ppf_out)
        optimizer.build_nvt(os.path.abspath(ppf_out))

    if cmd == 'optimize-npt':
        ppf_file = os.path.abspath(sys.argv[2])
        optimizer.optimize_npt(ppf_file)

    if cmd == 'npt':
        ppf_file = os.path.abspath(sys.argv[2])
        optimizer.run_npt(ppf_file)

    if cmd == 'plot':
        ppfs = sys.argv[2:]
        optimizer.plot(ppfs)
