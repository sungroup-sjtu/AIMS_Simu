import json
from flask import render_template, abort, request
from . import main
from .actions import StatAction
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


@main.route('/stat/', methods=['GET', 'POST'])
def show_stat():
    if request.method == 'GET':
        return render_template('stat.html')

    smiles_str = request.form.get('smiles')
    smiles_list = [line.strip() for line in smiles_str.splitlines()]
    smiles_list = list(filter(lambda x: x != '', smiles_list))

    statAction = StatAction()
    statAction.update_mol_task_list(smiles_list)

    ### NIST
    dens_Tm25_nist_list = statAction.get_density_nist(T='Tm25')
    dens_Tvap_nist_list = statAction.get_density_nist(T='Tvap')
    dens_Tcx8_nist_list = statAction.get_density_nist(T='Tcx8')

    hvap_Tm25_list = statAction.get_hvap_nist(T='Tm25')
    hvap_Tvap_list = statAction.get_hvap_nist(T='Tvap')
    hvap_Tcx8_list = statAction.get_hvap_nist(T='Tcx8')

    cp_Tm25_nist_list = statAction.get_cp_nist(T='Tm25')
    cp_Tvap_nist_list = statAction.get_cp_nist(T='Tvap')
    cp_Tcx8_nist_list = statAction.get_cp_nist(T='Tcx8')

    tc_nist_list = statAction.get_tc_nist()
    dc_nist_list = statAction.get_dc_nist()

    dgas_Tcx8_nist_list = statAction.get_dgas_nist(T='Tcx8')

    st_Tm25_nist_list = statAction.get_st_nist(T='Tm25')
    st_Tvap_nist_list = statAction.get_st_nist(T='Tvap')
    st_Tcx8_nist_list = statAction.get_st_nist(T='Tcx8')

    # ei_Tm25_nist_list = statAction.get_ei_nist(T='Tm25')
    # ei_Tvap_nist_list = statAction.get_ei_nist(T='Tvap')
    # ei_Tcx8_nist_list = statAction.get_ei_nist(T='Tcx8')

    # pv_Tm25_nist_list = statAction.get_pv_nist(T='Tm25')
    # pv_Tvap_nist_list = statAction.get_pv_nist(T='Tvap')
    # pv_Tcx8_nist_list = statAction.get_pv_nist(T='Tcx8')

    # sound_Tm25_nist_list = statAction.get_sound_nist(T='Tm25')
    # sound_Tvap_nist_list = statAction.get_sound_nist(T='Tvap')
    # sound_Tcx8_nist_list = statAction.get_sound_nist(T='Tcx8')

    ### Yaws
    # expan_Tm25_list = statAction.get_expansion(T='Tm25')
    # expan_Tvap_list = statAction.get_expansion(T='Tvap')
    # expan_Tcx8_list = statAction.get_expansion(T='Tcx8')

    import pylab
    import base64
    import numpy as np
    from io import BytesIO

    pngs = []
    for data_exp_sim_list in [
        dens_Tm25_nist_list, dens_Tvap_nist_list, dens_Tcx8_nist_list,
        hvap_Tm25_list, hvap_Tvap_list, hvap_Tcx8_list,
        cp_Tm25_nist_list, cp_Tvap_nist_list, cp_Tcx8_nist_list,
        tc_nist_list, dc_nist_list,
        dgas_Tcx8_nist_list,
        st_Tm25_nist_list, st_Tvap_nist_list, st_Tcx8_nist_list,
        # ei_Tm25_nist_list, ei_Tvap_nist_list, ei_Tcx8_nist_list,
        # pv_Tm25_nist_list, pv_Tvap_nist_list, pv_Tcx8_nist_list,
        # sound_Tm25_nist_list, sound_Tvap_nist_list, sound_Tcx8_nist_list,
        # expan_Tm25_list, expan_Tvap_list, expan_Tcx8_list,
    ]:
        ref_list = []
        sim_list = []
        dev_list = []
        u_list = []
        for ref, sim, u in data_exp_sim_list:
            abs_dev = abs(sim / ref - 1) * 100
            if abs_dev < 50:
                dev_list.append(abs_dev)  # incorrect expt. data
            ref_list.append(ref)
            sim_list.append(sim)
            u_list.append(u / ref * 100)

        p = BytesIO()
        fig = pylab.figure(figsize=(13, 4))
        fig.tight_layout()

        sp1 = fig.add_subplot(131)
        sp1.set_xlabel('Expt.')
        sp1.set_ylabel('Simu.')
        sp1.plot(ref_list, ref_list, '-')
        sp1.plot(ref_list, sim_list, '.', alpha=0.7, label='%i molecules' % len(data_exp_sim_list))
        sp1.legend()

        y, _ = np.histogram(dev_list, bins=30, range=[0, 30])
        x = (_[1:] + _[:-1]) / 2
        sp2 = fig.add_subplot(132)
        sp2.set_xlabel('Absolute deviation (%)')
        sp2.set_ylabel('Number of molecules')
        sp2.bar(x, y, color='C1', alpha=0.7, label='Mean = %.1f %%' % np.mean(dev_list))
        sp2.legend()

        y, _ = np.histogram(u_list, bins=30, range=[0, 30])
        x = (_[1:] + _[:-1]) / 2
        sp3 = fig.add_subplot(133)
        sp3.set_xlabel('Uncertainty of expt. data(%)')
        sp3.set_ylabel('Number of molecules')
        sp3.bar(x, y, alpha=0.7, label='Mean = %.1f %%' % np.mean(u_list))
        sp3.legend()

        pylab.savefig(p, format='png')
        pylab.close()
        pngs.append(base64.b64encode(p.getvalue()).decode())

    return render_template('stat.html', smiles_list=smiles_list, pngs=pngs)
