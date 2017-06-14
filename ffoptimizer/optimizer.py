#!/usr/bin/env python3
# coding=utf-8

import sys
import os
import time
import json
from collections import OrderedDict

import numpy as np
from lmfit import Parameters, Minimizer

sys.path.append('..')

from ffoptimizer.ppf import PPF
from ffoptimizer.db import DB
from ffoptimizer.models import Target, Result

from sqlalchemy import and_

CWD = os.getcwd()


class Optimizer():
    def __init__(self, db_file):
        self.db = DB(db_file)
        self.db.conn()

    @property
    def iteration(self):
        target = self.db.session.query(Target).first()
        return target.iteration

    def optimize_npt(self, ppf_file):
        LOG = os.path.join(CWD, 'Opt-%s.log' % os.path.basename(ppf_file)[:-4])

        def residual(params: Parameters):
            ### if result exist in database, ignore calculation
            result = self.db.session.query(Result).filter(
                and_(
                    Result.ppf == ppf_file,
                    Result.parameter == str(params)
                )).first()
            if result is not None:
                R = result.residual
                if R is not None:
                    return json.loads(R)
            ###

            ### save ppf file and run NPT
            ppf = PPF(ppf_file)
            paras = OrderedDict()
            for k, v in params.items():
                paras[k] = v.value
            ppf.set_lj_para(paras)
            ppf.write(ppf_file)

            for target in self.db.session.query(Target).all():
                if not target.npt_started(ppf_file):
                    ### next iteration
                    target.iteration += 1
                    target.run_npt(ppf_file, paras)
            self.db.session.commit()
            ###

            while True:
                FINISHED = True
                for target in self.db.session.query(Target).all():
                    if not target.npt_finished(ppf_file):
                        FINISHED = False
                        break

                if FINISHED:
                    break
                else:
                    current_time = time.strftime('%m-%d %H:%M')
                    print(current_time + ' Job still running. Wait ...')
                    time.sleep(60)

            R = []
            targets = self.db.session.query(Target).all()
            for target in targets:
                dens, hvap = target.get_npt_result(os.path.basename(ppf_file)[:-4])
                R.append((dens - target.density) / target.density * 100 * target.wDensity)  # deviation  percent
                R.append((hvap - target.hvap) / target.hvap * 100 * target.wHvap)  # deviation percent

            ### save result to database
            if result is None:
                result = Result(ppf=ppf_file, parameter=str(params))

            result.residual = json.dumps(R)
            self.db.session.add(result)
            self.db.session.commit()
            ###

            ### write current parameters and residual to log
            txt = '\nITERATION: %i\n' % self.iteration
            txt += '\nPARAMETERS:\n'
            for k, v in params.items():
                txt += '%10.5f %s\n' % (v.value, k)
            txt += '\nRESIDUE:\n'
            for i, r in enumerate(R):
                name_prop = targets[i//2].name
                name_prop += ' density' if i % 2 == 0 else ' hvap'
                txt += '%8.2f %s\n' % (r, name_prop)
            txt += '\nRSQ: %.2f\n' % np.sum(list(map(lambda x: x ** 2, R)))

            print(txt)
            with open(LOG, 'a') as log:
                log.write(txt)
            ###

            return R

        def jacobian(params: Parameters):
            ### if result exist in database, ignore calculation
            result = self.db.session.query(Result).filter(
                and_(
                    Result.ppf == ppf_file,
                    Result.parameter == str(params)
                )).first()
            if result is not None:
                J = result.jacobian
                if J is not None:
                    return json.loads(J)
            ###

            paras = OrderedDict()
            for k, v in params.items():
                paras[k] = v.value

            J = []
            targets = self.db.session.query(Target).all()
            for target in targets:
                dDens, dHvap = target.get_dDens_dHvap_from_paras(ppf_file, paras)
                J.append([i / target.density * 100 * target.wDensity for i in dDens])  # deviation  percent
                J.append([i / target.hvap * 100 * target.wHvap for i in dHvap])  # deviation  percent

            ### save result to database
            if result is None:
                result = Result(ppf=ppf_file, parameter=str(params))

            result.jacobian = json.dumps(J)
            self.db.session.add(result)
            self.db.session.commit()
            ###

            ### write Jacobian to log
            txt = '\nJACOBIAN MATRIX:\n'
            for k in params.keys():
                txt += '%9s' % k
            txt += '\n'
            for i, row in enumerate(J):
                name_prop = targets[i//2].name
                name_prop += ' density' if i % 2 == 0 else ' hvap'
                for item in row:
                    txt += '%9.2f' % item
                txt += ' %s\n' % name_prop

            print(txt)
            with open(LOG, 'a') as log:
                log.write(txt)
            ###

            return J

        def callback(params: Parameters, iter: int, res: [float]):
            ### clear _finished_ and job.sh for next iteration
            for target in self.db.session.query(Target).all():
                target.clear_npt_result(ppf_file)

        ppf = PPF(ppf_file)
        params = Parameters()
        for k, v in ppf.adj_lj_paras.items():
            if k.endswith('r0'):
                params.add(k, value=v, min=2, max=5)
            elif k.endswith('e0'):
                params.add(k, value=v, min=0.005, max=0.2)
            elif k.endswith('bi'):
                params.add(k, value=v, min=-0.3, max=0)

        minimize = Minimizer(residual, params, iter_cb=callback)
        result = minimize.leastsq(Dfun=jacobian, ftol=0.0001)
        print(result.lmdif_message, '\n')

        return result.params

    def npt(self, ppf_file):
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

    elif cmd == 'optimize':
        ppf_file = os.path.abspath(sys.argv[2])
        params_out = optimizer.optimize_npt(ppf_file)

    elif cmd == 'npt':
        ppf_file = os.path.abspath(sys.argv[2])
        optimizer.npt(ppf_file)

    elif cmd == 'plot':
        ppfs = sys.argv[2:]
        optimizer.plot(ppfs)

    else:
        print('Unknown command')
