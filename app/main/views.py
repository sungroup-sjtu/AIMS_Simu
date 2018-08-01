import json
from flask import render_template, abort, request
from . import main
from .actions import StatAction, VleAction
from ..models import Compute, Task, Job


@main.route('/')
def index():
    return render_template('index.html')


@main.route('/compute/<int:compute_id>')
def show_compute(compute_id):
    compute = Compute.query.get(compute_id)
    return render_template('compute.html', compute=compute)


@main.route('/task/<int:task_id>')
def show_task(task_id):
    task = Task.query.get(task_id)
    return render_template('task.html', task=task)


@main.route('/job/<int:job_id>')
def show_job(job_id):
    job = Job.query.get(job_id)
    return render_template('job.html', job=job)


def get_pngs_from_data(k_data_exp_sim_list):
    import pylab
    import base64
    import numpy as np
    from io import BytesIO

    pngs = []
    for k, data_exp_sim_list in k_data_exp_sim_list:
        ref_list = []
        sim_list = []
        dev_list = []
        absdev_list = []
        u_list = []
        for ref, sim, u in data_exp_sim_list:
            dev = (sim / ref - 1) * 100
            if abs(dev) < 50:  # incorrect expt. data
                dev_list.append(dev)
                absdev_list.append(abs(dev))
                ref_list.append(ref)
                sim_list.append(sim)
                u_list.append(u / ref * 100)

        p = BytesIO()
        fig = pylab.figure(figsize=(13, 4))

        sp1 = fig.add_subplot(131)
        sp1.set_xlabel('Expt.')
        sp1.set_ylabel('Simu.')
        sp1.plot(ref_list, ref_list, '-')
        sp1.plot(ref_list, sim_list, '.', alpha=0.7)

        props = dict(boxstyle='round', facecolor='wheat', alpha=0.5)
        text = '%s\n%i molecules\nMDEV = %.1f %%\nSTD = %.1f %%\nMUD = %.1f %%' \
               % (k, len(ref_list), np.mean(dev_list), np.std(dev_list), np.mean(absdev_list))
        sp1.text(0.05, 0.95, text, transform=sp1.transAxes, va='top', bbox=props)

        y, _ = np.histogram(absdev_list, bins=30, range=[0, 30])
        x = (_[1:] + _[:-1]) / 2
        sp2 = fig.add_subplot(132)
        sp2.set_xlabel('Unsigned deviation (%)')
        sp2.set_ylabel('Number of molecules')
        sp2.bar(x, y, color='C1', alpha=0.7)

        text = 'Unsigned deviation to expt. data\n%s\n%i molecules\nMUD = %.1f %%' \
               % (k, len(ref_list), np.mean(absdev_list))
        sp2.text(0.95, 0.95, text, transform=sp2.transAxes, va='top', ha='right', bbox=props)

        y, _ = np.histogram(u_list, bins=30, range=[0, 30])
        x = (_[1:] + _[:-1]) / 2
        sp3 = fig.add_subplot(133)
        sp3.set_xlabel('Uncertainty (%)')
        sp3.set_ylabel('Number of molecules')
        sp3.bar(x, y, alpha=0.7)

        text = 'Uncertainty of expt. data\n%s\n%i molecules\nMean uncertainty = %.1f %%' \
               % (k, len(ref_list), np.mean(u_list))
        sp3.text(0.95, 0.95, text, transform=sp3.transAxes, va='top', ha='right', bbox=props)

        pylab.savefig(p, format='png')
        pylab.close()
        pngs.append(base64.b64encode(p.getvalue()).decode())

    return pngs


@main.route('/stat/pvt', methods=['GET', 'POST'])
def show_stat_pvt():
    if request.method == 'GET':
        return render_template('stat.html')

    smiles_str = request.form.get('smiles')
    smiles_list = [line.strip() for line in smiles_str.splitlines()]
    smiles_list = list(set(filter(lambda x: x != '', smiles_list)))
    print(len(smiles_list))

    statAction = StatAction()
    statAction.update_mol_task_list(smiles_list)
    print(len(statAction.nist_list))

    dens_Tm25_nist_list = statAction.get_density_nist(_T='Tm25')
    dens_Tvap_nist_list = statAction.get_density_nist(_T='Tvap')
    dens_Tcx8_nist_list = statAction.get_density_nist(_T='Tcx8')
    hvap_Tm25_list = statAction.get_hvap_nist(_T='Tm25')
    hvap_Tvap_list = statAction.get_hvap_nist(_T='Tvap')
    hvap_Tcx8_list = statAction.get_hvap_nist(_T='Tcx8')
    cp_Tm25_nist_list = statAction.get_cp_nist(_T='Tm25')
    cp_Tvap_nist_list = statAction.get_cp_nist(_T='Tvap')
    cp_Tcx8_nist_list = statAction.get_cp_nist(_T='Tcx8')
    # sound_Tm25_nist_list = statAction.get_sound_nist(T='Tm25')
    # sound_Tvap_nist_list = statAction.get_sound_nist(T='Tvap')
    # sound_Tcx8_nist_list = statAction.get_sound_nist(T='Tcx8')

    k_data_exp_sim_list = [
        ('density @ Tm+', dens_Tm25_nist_list),
        ('density @ Tvap', dens_Tvap_nist_list),
        ('density @ T0.8*', dens_Tcx8_nist_list),
        ('Hvap @ Tm+', hvap_Tm25_list),
        ('Hvap @ Tvap', hvap_Tvap_list),
        ('Hvap @ T0.8*', hvap_Tcx8_list),
        ('Cp @ Tm+', cp_Tm25_nist_list),
        ('Cp @ Tvap', cp_Tvap_nist_list),
        ('Cp @ T0.8*', cp_Tcx8_nist_list),
        # sound_Tm25_nist_list, sound_Tvap_nist_list, sound_Tcx8_nist_list,
    ]
    pngs = get_pngs_from_data(k_data_exp_sim_list)

    return render_template('stat.html', smiles_list=smiles_list, pngs=pngs)


@main.route('/stat/vle', methods=['GET', 'POST'])
def show_stat_vle():
    if request.method == 'GET':
        return render_template('stat.html')

    smiles_str = request.form.get('smiles')
    smiles_list = [line.strip() for line in smiles_str.splitlines()]
    smiles_list = list(set(filter(lambda x: x != '', smiles_list)))
    print(len(smiles_list))

    statAction = StatAction()
    statAction.update_mol_slab_list(smiles_list)
    print(len(statAction.nist_list))

    tc_nist_list = statAction.get_tc_nist()
    dc_nist_list = statAction.get_dc_nist()
    dgas_Tcx8_nist_list = statAction.get_dgas_nist(_T='Tcx8')
    st_Tm25_nist_list = statAction.get_st_nist(_T='Tm25')
    st_Tvap_nist_list = statAction.get_st_nist(_T='Tvap')
    st_Tcx8_nist_list = statAction.get_st_nist(_T='Tcx8')

    k_data_exp_sim_list = [
        ('critical temperatrue', tc_nist_list),
        ('critical density', dc_nist_list),
        ('density of vapor @ T0.8*', dgas_Tcx8_nist_list),
        ('surface tension @ Tm+', st_Tm25_nist_list),
        ('surface tension @ Tvap', st_Tvap_nist_list),
        ('surface tension @ T0.8*', st_Tcx8_nist_list),
    ]
    pngs = get_pngs_from_data(k_data_exp_sim_list)

    return render_template('stat.html', smiles_list=smiles_list, pngs=pngs)


@main.route('/vle/', methods=['GET', 'POST'])
def show_vle():
    if request.method == 'GET':
        return render_template('vle.html')

    smiles_str = request.form.get('smiles')
    smiles_list = [line.strip() for line in smiles_str.splitlines()]
    smiles_list = list(filter(lambda x: x != '', smiles_list))

    import pylab
    import base64
    from io import BytesIO
    pngs = []

    dg_absdev_list = []
    dg_dev_list = []
    dl_absdev_list = []
    dl_dev_list = []
    dl_u_list = []
    dg_u_list = []

    action = VleAction()
    for smiles in smiles_list:
        if not action.get_vle(smiles):
            continue

        dl = action.dl_list[-1]
        dl_u = action.dl_u_list[-1]
        dls = action.dls_list[-1]
        dg = action.dg_list[-1]
        dg_u = action.dg_u_list[-1]
        dgs = action.dgs_list[-1]
        if dl != None and dls != None:
            dl_absdev_list.append(abs(dls - dl))
            dl_dev_list.append(abs(dls - dl) / dl * 100)
            dl_u_list.append(dl_u / dl * 100)
        if dg != None and dgs != None:
            dg_absdev_list.append(abs(dgs - dg))
            dg_dev_list.append(abs(dgs - dg) / dg * 100)
            dg_u_list.append(dg_u / dg * 100)

        p = BytesIO()
        fig = pylab.figure(figsize=(2.7, 2.7))
        fig.tight_layout()

        sp1 = fig.add_subplot(111)
        sp1.plot(action.dl_list, action.t_list, '--', linewidth=2, color='C0')
        sp1.plot(action.dg_list, action.t_list, '--', linewidth=2, color='C0')
        sp1.plot(action.dls_list, action.ts_list, '.', markersize=10, alpha=0.8, color='C1')
        sp1.plot(action.dgs_list, action.ts_list, '.', markersize=10, alpha=0.8, color='C1')

        props = dict(boxstyle='round', facecolor='wheat', alpha=0.5)
        text = '%s\n%s' % (action.nist.formula, action.nist.name)
        sp1.text(0.5, 0.05, text, transform=sp1.transAxes, ha='center', va='bottom', bbox=props)

        pylab.savefig(p, format='png')
        pylab.close()
        pngs.append(base64.b64encode(p.getvalue()).decode())

    import numpy as np
    print(len(dl_dev_list), len(dg_dev_list))
    print(np.mean(dl_absdev_list), np.mean(dl_dev_list), np.mean(dl_u_list))
    print(np.mean(dg_absdev_list), np.mean(dg_dev_list), np.mean(dg_u_list))

    return render_template('vle.html', smiles_list=smiles_list, pngs=pngs)
