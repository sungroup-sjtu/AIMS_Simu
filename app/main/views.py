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


@main.route('/stat/')
def show_stat():
    return render_template('stat.html')


@main.route('/stat/', methods=['POST'])
def get_stat_json():
    smiles_str = request.form.get('smiles')
    smiles_list = [line.strip() for line in smiles_str.splitlines()]
    smiles_list = list(filter(lambda x: x != '', smiles_list))

    statAction = StatAction()
    statAction.update_mol_task_list(smiles_list)

    dens_298K_list = statAction.get_density(smiles_list, T=298)
    dens_Tvap_list = statAction.get_density(smiles_list, T='Tvap')
    dens_Tc_list = statAction.get_density(smiles_list, T='Tc')

    Cp_298K_list = statAction.get_Cp(smiles_list, T=298)
    Cp_Tvap_list = statAction.get_Cp(smiles_list, T='Tvap')
    Cp_Tc_list = statAction.get_Cp(smiles_list, T='Tc')

    Hvap_298K_list = statAction.get_Hvap(smiles_list, T=298)
    Hvap_Tvap_list = statAction.get_Hvap(smiles_list, T='Tvap')

    return json.dumps(dict(n_smiles=len(smiles_list),
                           n_matched=len(statAction.yaws_list),
                           datasets=[
                               dict(name='dens_298K', data=dens_298K_list),
                               dict(name='dens_Tvap', data=dens_Tvap_list),
                               dict(name='dens_Tc', data=dens_Tc_list),
                               dict(name='Cp_298K', data=Cp_298K_list),
                               dict(name='Cp_Tvap', data=Cp_Tvap_list),
                               dict(name='Cp_Tc', data=Cp_Tc_list),
                               dict(name='Hvap_298K', data=Hvap_298K_list),
                               dict(name='Hvap_Tvap', data=Hvap_Tvap_list),
                           ]
                           )
                      )
