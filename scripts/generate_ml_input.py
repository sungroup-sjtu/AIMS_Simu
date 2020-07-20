#!/usr/bin/env python3
# coding=utf-8

import sys
import pandas as pd

sys.path.append('..')
from app import create_app
from app.models import *
from app.models_nist import *
from app.models_cv import Cv
from app.models_ilthermo import *
from app.selection import *
from mstools.analyzer.fitting import (
    VTFfit,
    fit_vle_dminus,
    fit_vle_dplus,
    fit_vle_pvap,
    fit_vle_st
)


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


def main():
    import argparse
    parser = argparse.ArgumentParser(description='This is a code to generate input txt from experimental data')
    parser.add_argument('-t', '--type', type=str, help='The type of the output(SIM, EXP)')
    parser.add_argument('--database', type=str, help='The type of the experimental database(ILTHERMO, NIST)')
    parser.add_argument('--procedure', type=str, help='The procedure of simulation(npt)')
    parser.add_argument('--selection', type=bool, help='Use the selection function', default=True)
    parser.add_argument('--nheavy', type=int, help='Only 5 < heavy atoms < nheavy is considered', default=16)
    parser.add_argument('--training', type=str, help='The txt file for training sets smiles', default=None)
    parser.add_argument('--errormolecules', type=bool,
                        help='Find the molecules in training set without simulation data',
                        default=False)
    parser.add_argument('--rawdata', type=bool, help='Use the raw simulation data, otherwise fit data', default=False)
    parser.add_argument('--explicit', type=bool, help='Only use the simulation data with unambiguous smiles',
                        default=True)

    args = parser.parse_args()
    app = create_app(args.procedure)
    app.app_context().push()
    if args.type == 'EXP' and args.database == 'NIST':
        def get_exp_data(property, uncertainty=True, T=False):
            df = pd.DataFrame({'inchi': [], 'SMILES': [], 'nheavy': []})
            if T:
                df['T'] = []
                df['rel_T'] = []
                df['tt'] = []
                df['tc'] = []
            df[property] = []
            if uncertainty:
                df['%s_u' % property] = []

            molecules = NistMolecule.query.filter(NistMolecule.n_heavy >
                                                  1).limit(5)
            # 5).filter(NistMolecule.n_heavy < 16)
            for i, mol in enumerate(molecules):
                if mol.remark != 'selected':
                    continue
                sys.stdout.write('\r%i / %i. %s\t\t\t\t' % (
                i, molecules.count(), mol.smiles))
                if T:
                    if mol.tt is None or mol.tc is None:
                        continue
                    datas_mol = mol.datas
                    datas = datas_mol.filter(
                        NistData.property == NistProperty.query.filter(
                            NistProperty.name == property).first())
                    if datas is not None:
                        for data in datas:
                            rel_T = (data.t - mol.tt) / (mol.tc - mol.tt)
                            df.loc[df.shape[0]] = mol.inchi, mol.smiles, \
                                                  mol.n_heavy, data.t, rel_T, \
                                                  mol.tt, mol.tc, \
                                                  data.value, data.uncertainty
                else:
                    if property == 'tb':
                        v, u = mol.tb, mol.tb_u
                    elif property == 'tt':
                        v, u = mol.tt, mol.tt_u
                    elif property == 'tc':
                        v, u = mol.tc, mol.tc_u
                    elif property == 'pc':
                        v, u = mol.pc, mol.pc_u
                    elif property == 'dc':
                        v, u = mol.dc, mol.dc_u
                    elif property == 'hfus':
                        v, u = mol.hfus, mol.hfus_u
                    elif property == 'st':
                        spline = mol.splines.filter(
                            NistSpline.property_id == 8).first()
                        if mol.tb is not None and spline is not None:
                            v, u = spline.get_data(T=mol.tb)
                        else:
                            v, u = None, None
                    elif property == 'vis':
                        spline = mol.splines.filter(
                            NistSpline.property_id == 7).first()
                        if mol.tb is not None and spline is not None:
                            v, u = spline.get_data(T=mol.tb)
                        else:
                            v, u = None, None
                    else:
                        v, u = None, None
                    if v is not None and u is not None:
                        df.loc[df.shape[0]] = mol.inchi, mol.smiles, \
                                              mol.n_heavy, v, u
            df.to_csv('%s.txt' % property, sep=' ', index=False)

        def get_exp_data_fitcoef(property, repeat=100):
            df = pd.DataFrame({'inchi': [], 'SMILES': []})
            if property == 'viscosity-lg':
                df['c1'] = []
                df['c2'] = []
                df['c3'] = []
                df['c1_u'] = []
                df['c2_u'] = []
                df['c3_u'] = []
            molecules = NistMolecule.query.filter(
                NistMolecule.n_heavy > 5).filter(NistMolecule.n_heavy < 16)
            for i, mol in enumerate(molecules):
                if mol.remark != 'selected':
                    continue
                sys.stdout.write('\r%i / %i. %s\t\t\t\t' % (
                i, molecules.count(), mol.smiles))
                datas_mol = mol.datas
                datas = datas_mol.filter(
                    NistData.property == NistProperty.query.filter(
                        NistProperty.name == property).first())
                if datas.count() < 10:
                    continue
                t_list = []
                v_list = []
                u_list = []
                if datas is not None:
                    for data in datas:
                        if data.value is None or data.uncertainty is None:
                            continue
                        if data.uncertainty / data.value > 1.0:
                            continue
                        t_list.append(data.t)
                        v_list.append(data.value)
                        u_list.append(data.uncertainty)
                if property == 'viscosity-lg':
                    # print(len(v_list), len(u_list))
                    coef_list = []
                    j = 0
                    k = 0
                    while j < repeat:
                        k += 1
                        if k == repeat * 10:
                            break
                        if len(v_list) < 5:
                            k = repeat * 10
                            break
                        data = np.random.normal(
                            loc=v_list,
                            scale=u_list,
                            size=len(v_list)
                        )
                        if data.min() < 0:
                            continue
                        try:
                            coef, score = VTFfit(t_list, data)
                        except:
                            k = repeat * 10
                            break
                        else:
                            coef_list.append(coef)
                            j += 1
                    if k == repeat * 10:
                        continue
                    coef_v = np.mean(coef_list, axis=0)
                    coef_u = np.std(coef_list, axis=0)
                    df.loc[df.shape[0]] = mol.inchi, mol.smiles, coef_v[0], \
                                          coef_v[1], \
                                          coef_v[2], coef_u[0], coef_u[1], \
                                          coef_u[2],
            df.to_csv('%s.txt' % property, sep=' ', index=False)

        def get_vle_coefs(n=100):
            df = pd.DataFrame({'inchi': [], 'SMILES': [], 'tc': [], 'dc': [],
                               'A': [], 'B': [], 'tmin': [], 'tmax': [],
                               'score_minus': [], 'score_plus': []})
            splines_liq = NistSpline.query.filter(NistSpline.property_id == 2)
            splines_gas = NistSpline.query.filter(NistSpline.property_id == 3)
            molecules = NistMolecule.query.filter(NistMolecule.n_heavy > 1)
            for i, mol in enumerate(molecules):
                if mol.remark != 'selected':
                    continue
                sys.stdout.write('\r%i / %i. %s\t\t\t\t' % (
                i, molecules.count(), mol.smiles))
                sp_liq = splines_liq.filter(
                    NistSpline.molecule == mol
                ).first()
                if sp_liq is None:
                    continue
                sp_gas = splines_gas.filter(
                    NistSpline.molecule == mol
                ).first()
                if sp_gas is None:
                    continue
                t_min = max(sp_liq.t_min, sp_gas.t_min)
                t_max = min(sp_liq.t_max, sp_gas.t_max)
                t_list = np.linspace(t_min, t_max, n)
                d_liq = sp_liq.get_data(t_list)[0]
                d_gas = sp_gas.get_data(t_list)[0]
                try:
                    d_minus = d_liq - d_gas
                    coef, score_minus = fit_vle_dminus(t_list, d_minus)
                    tc = coef[0]
                    B = coef[1]
                    d_plus = d_liq + d_gas
                    coef, score_plus = fit_vle_dplus(t_list, d_plus, tc)
                    dc = coef[0]
                    A = coef[1]
                    df.loc[df.shape[0]] = mol.inchi, mol.smiles, tc, dc, A, B, \
                                          t_min, t_max, score_minus, score_plus
                except:
                    continue
            df.to_csv('vle-coef.txt', sep=' ', index=False)

        def get_coefs(property_id):
            df = pd.DataFrame({'inchi': [], 'SMILES': [], 'nheavy': [],
                               'score': []})
            if property_id == 1:
                df['c1'] = []
                df['c2'] = []
            elif property_id == 7:
                df['c1'] = []
                df['c2'] = []
                df['c3'] = []
            elif property_id == 8:
                df['A'] = []
                df['n'] = []
                df['tc'] = []
            molecules = NistMolecule.query.filter(NistMolecule.n_heavy > 1)
            for i, mol in enumerate(molecules):
                if mol.remark != 'selected':
                    continue
                sys.stdout.write('\r%i / %i. %s\t\t\t\t' % (
                i, molecules.count(), mol.smiles))
                datas = mol.datas.filter(NistData.property_id == property_id)
                if datas.count() < 10:
                    continue
                t_list = []
                v_list = []
                u_list = []
                for data in datas:
                    if data.value is None or data.uncertainty is None:
                        continue
                    t_list.append(data.t)
                    v_list.append(data.value)
                    u_list.append(data.uncertainty)
                try:
                    if property_id == 1:
                        coef, score = fit_vle_pvap(t_list, v_list)
                        df.loc[df.shape[0]] = mol.inchi, mol.smiles, \
                                              mol.n_heavy, score, coef[0], \
                                              coef[1]
                    elif property_id == 7:
                        coef, score = VTFfit(t_list, v_list)
                        df.loc[df.shape[0]] = mol.inchi, mol.smiles, \
                                              mol.n_heavy, score, coef[0], \
                                              coef[1], coef[2]
                    elif property_id == 8 and mol.tc is not None:
                        coef, score = fit_vle_st(t_list, v_list, mol.tc)
                        df.loc[df.shape[0]] = mol.inchi, mol.smiles, \
                                              mol.n_heavy, score, coef[0], \
                                              coef[1], mol.tc
                except:
                    continue
            df.to_csv('coef-%i.txt' % property_id, sep=' ', index=False)
        # get_exp_data('tt')
        #get_exp_data('tb')
        #get_exp_data('tc')
        #get_exp_data('pc')
        #get_exp_data('dc')
        #get_exp_data('hfus')
        get_exp_data('pvap-lg', T=True)
        get_exp_data('density-lg', T=True)
        get_exp_data('density-gl', T=True)
        # get_exp_data('hvap-lg', T=True)
        # get_exp_data('cp-lg', T=True)
        # get_exp_data('sound-lg', T=True)
        #get_exp_data('viscosity-lg', T=True)
        # get_exp_data('st-lg', T=True)
        #get_coefs(1)
        #get_coefs(7)
        #get_coefs(8)
        # get_exp_data_fitcoef('viscosity-lg')
        # get_exp_data_fitcoef_old('viscosity-lg')
        # get_vle_coefs(100)
    elif args.type == 'EXP' and args.database == 'ILTHERMO':
        def get_exp_vis_coef():
            df = pd.DataFrame({'SMILES': [], 'c1': [], 'c2': [], 'c3': [],
                               'tmin': [], 'tmax': []})
            splines = Spline.query.filter(Spline.property_id == 50)
            for spline in splines:
                smiles = spline.molecule.smiles()
                if smiles is None:
                    continue
                coefs, score = json.loads(spline.coef_VTF)
                df.loc[df.shape[0]] = smiles, coefs[0], coefs[1], coefs[2], \
                                      spline.t_min, spline.t_max
            df.to_csv('vis-coefs.txt', sep=' ', index=False)

        def get_exp_data(p_id):
            df = pd.DataFrame({'SMILES': [], 'T': [], 'P': [], 'vis': [],
                               'vis_u': []})
            datas = Data.query.filter(Data.property_id == p_id)
            for data in datas:
                smiles = data.molecule.smiles()
                p = 100 if data.p is None else data.p
                if smiles is None:
                    continue
                df.loc[df.shape[0]] = smiles, data.t, p, data.value, data.stderr
            df.to_csv('vis.txt', sep=' ', index=False)

        get_exp_vis_coef()
        get_exp_data(50)
    elif args.type == 'SIM' and args.errormolecules and args.training:
        info = pd.read_csv(args.training, sep='\s+', header=0)
        training_smiles_list = []
        for smiles in info['SMILES'].unique():
            training_smiles_list.append(get_canonical_smiles(smiles))
        print('%i molecules in training set' % len(training_smiles_list))
        for smiles in training_smiles_list:
            task = Task.query.filter(Task.procedure == args.procedure).filter(
                Task.smiles_list == json.dumps([smiles])).first()
            if task is None:
                print('There is no task for %s' % smiles)
                continue
            if task.status in [Compute.Status.FAILED, Compute.Status.STARTED]:
                print('Task for %s not finished' % smiles)
                continue
            if not task_selection(task, select=args.selection):
                print('%s does not satisfy selection rule' % smiles)
                continue
            post_result = task.get_post_result()
            if post_result is None:
                print('There is no post_result for %s' % smiles)
                continue
            t = task.get_t_list()[0]
            post_data = task.get_post_data(T=t)
            if args.procedure == 'npt':
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
            elif args.procedure == 'nvt-slab':
                st = post_data.get('st')
                tc = post_data.get('tc')
                dc = post_data.get('dc')
                if st is None:
                    print('%s has no st' % json.loads(task.smiles_list))
                if tc is None:
                    print('%s has no tc' % json.loads(task.smiles_list))
                if dc is None:
                    print('%s has no dc' % json.loads(task.smiles_list))
    elif args.type == 'SIM' and args.procedure == 'npt':
        if args.training is not None:
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
            info = pd.read_csv(args.training, sep='\s+', header=0)
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
        tasks = Task.query.filter(Task.procedure == args.procedure)
        for task in tasks:
            if task.status in [Compute.Status.FAILED, Compute.Status.STARTED]:
                continue
            if not task_selection(task, select=args.selection):
                continue
            if args.explicit and has_stereo_isomer(json.loads(task.smiles_list)[0]):
                continue
            post_result = task.get_post_result()
            if post_result is None or post_result['density-poly4'][-1] < 0.999:
                continue
            smiles = task.get_smiles_list()[0]
            py_mol = task.get_mol_list()[0]
            if not 5 < get_heavy_atom_numbers(smiles) < args.nheavy:
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
                    if args.rawdata and job_result.get('density') is not None:
                        density = job_result.get('density')[0]
                    else:
                        density = post_data.get('density')
                    if args.rawdata and job_result.get('einter') is not None:
                        einter = job_result.get('einter')[0] / min(n_mol_list)
                    else:
                        einter = post_data.get('einter')
                    if args.rawdata and job_result.get('expansion') is not None:
                        expansion = job_result.get('expansion')[0]
                    else:
                        expansion = post_data.get('expansion')
                    if args.rawdata and job_result.get('compress') is not None:
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
    elif args.type == 'SIM' and args.procedure == 'nvt-slab':
        if args.training is not None:
            f_train1 = open('result-ML-st-train.txt', 'w')
            f_train2 = open('result-ML-tc-train.txt', 'w')
            f_train3 = open('result-ML-dc-train.txt', 'w')
            f_train1.write('SMILES T st\n')
            f_train2.write('SMILES tc\n')
            f_train3.write('SMILES dc\n')
            info = pd.read_csv(args.training, sep='\s+', header=0)
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
        tasks = Task.query.filter(Task.procedure == args.procedure)
        for task in tasks:
            if task.status in [Compute.Status.FAILED, Compute.Status.STARTED]:
                continue
            if not task_selection(task, select=args.selection):
                continue
            if args.explicit and has_stereo_isomer(json.loads(task.smiles_list)[0]):
                continue
            post_result = task.get_post_result()
            if post_result is None:
                continue
            smiles = task.get_smiles_list()[0]
            if not 5 < get_heavy_atom_numbers(smiles) < args.nheavy:
                continue
            jobs = task.jobs
            for i, t in enumerate(task.get_t_list()):
                job_result = json.loads(jobs.filter(Job.t == t).first().result)
                post_data = task.get_post_data(T=t)
                if args.rawdata and job_result.get('st') is not None:
                    st = job_result.get('st')[0]
                else:
                    st = post_data.get('st')

                write_data(smiles, training_smiles_list, t, None, st, f1, f_train1, positive=True)

                if i == 0:
                    tc = post_data.get('tc')
                    write_data(smiles, training_smiles_list, None, None, tc, f2, f_train2, positive=True)
                    dc = post_data.get('dc')
                    write_data(smiles, training_smiles_list, None, None, dc, f3, f_train3, positive=True)


if __name__ == '__main__':
    main()
