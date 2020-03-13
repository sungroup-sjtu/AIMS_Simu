import json
import re
import shutil
import sys
import time
import traceback
from datetime import datetime
from functools import partial

from sqlalchemy import Column, ForeignKey, Integer, Text, String, Boolean, DateTime, and_

from flask import current_app

from config import Config
from . import db

sys.path.append(Config.MS_TOOLS_DIR)
from mstools.simulation.procedure import Procedure
from mstools.wrapper import GMX
from mstools.utils import *
from mstools.smiles.smiles import *
import numpy as np, pandas as pd


NotNullColumn = partial(Column, nullable=False)


def init_simulation(procedure, extend=False, bugfix=False):
    from mstools.simulation import gmx as simulationEngine
    kwargs = {'packmol'   : current_app.packmol,
              'dff'       : current_app.dff,
              'gmx'       : current_app.gmx,
              'jobmanager': current_app.jobmanager
              }
    if extend:
        kwargs.update({
            'gmx'       : current_app.gmx_extend,
            'jobmanager': current_app.jm_extend
        })
    if bugfix:
        kwargs.update({
            'gmx_bin'   : current_app.config['BUGFIX_GMX_BIN'],
            'gmx_mdrun' : current_app.config['BUGFIX_GMX_MDRUN'],
            'jobmanager': current_app.jm_bugfix
        })

    if procedure in ['npt', 'npt-multi', 'npt-v-rescale', 'npt-2', 'npt-3']:
        sim = simulationEngine.Npt(**kwargs)
    elif procedure in ['nvt', 'nvt-multi', 'nvt-multi-2', 'nvt-multi-3']:
        sim = simulationEngine.Nvt(**kwargs)
    elif procedure == 'nvt-slab':
        sim = simulationEngine.NvtSlab(**kwargs)
    elif procedure == 'ppm':
        sim = simulationEngine.NptPPM(**kwargs)
    else:
        raise Exception('Unknown simulation procedure')

    sim.gmx._DIELECTRIC = Config.CHARGE_SCALE ** -2
    return sim


def wrapper_cls_func(cls_func_args_kwargs):
    cls, func, args, kwargs = cls_func_args_kwargs
    func = getattr(cls, func)
    return func(*args, **kwargs)


class PbsJob(db.Model):
    '''
    A PbsJob record corresponds to a PBS Job submitted to queue system like Slurm and Torque
    '''
    __tablename__ = 'pbs_job'
    id = NotNullColumn(Integer, primary_key=True)
    name = NotNullColumn(String(200))
    sh_file = NotNullColumn(Text)
    extend = NotNullColumn(Boolean, default=False)
    submitted = NotNullColumn(Boolean, default=False)
    bugfix = NotNullColumn(Boolean, default=False)

    def __repr__(self):
        return '<PbsJob: %i: %s %s>' % (self.id, self.name, self.submitted)

    @property
    def jm(self):
        if self.extend:
            return current_app.jm_extend
        elif self.bugfix:
            return current_app.jm_bugfix
        else:
            return current_app.jobmanager

    def submit(self, **kwargs):
        if self.jm.submit(self.sh_file, **kwargs):
            self.submitted = True
        else:
            self.submitted = False
            current_app.logger.warning('Submit PBS job failed %s' % self.sh_file)
        db.session.commit()

    def kill(self):
        if not self.jm.kill_job(self.name):
            current_app.logger.warning('Kill PBS job failed %s' % self.name)

    def is_running(self):
        '''
        If the pbs job missing, return False
        If the pbs job has finished, return False
        Otherwise, return True
        :return:
        '''
        return self.jm.is_running(self.name)


class Compute(db.Model):
    '''
    A Compute record corresponds to a high throughput computation submitted
    It contains a lot of molecules and physical states
    A Compute record can produce a lot of Task records
    '''
    __tablename__ = 'compute'
    id = NotNullColumn(Integer, primary_key=True)
    web_id = NotNullColumn(Integer)
    web_user_id = NotNullColumn(Integer)
    web_ip = NotNullColumn(String(200))
    time = NotNullColumn(DateTime, default=datetime.now)
    json_detail = NotNullColumn('json', Text)
    remark = NotNullColumn(Text)

    tasks = db.relationship('Task', lazy='dynamic')

    def __repr__(self):
        return '<Compute: %i: %s>' % (self.id, self.remark)

    def create_tasks(self, alter_task=False, unambiguity=True, create_task=True):
        # alter_task: alter the t_list and p_list of exist task
        # unambiguity: if the smiles is ambiguous (do not assigned stereo-conformation), do not create task
        current_app.logger.info('Create tasks from %s' % self)
        try:
            detail = json.loads(self.json_detail)
            procedure = detail['procedures'][0]
            combinations = detail['combinations']

            # check for procedure and prior
            if procedure not in Procedure.choices:
                raise Exception('Invalid simulation procedure: %s' % procedure)

            prior_procedure = Procedure.prior.get(procedure)

            N = 0
            M = 0
            for combination in combinations:
                # check for smiles and name
                smiles_list = combination['smiles']
                if unambiguity and True in [has_stereo_isomer(smiles) for smiles in smiles_list]:
                    continue

                name_list = combination.get('names')
                t_list = list(map(int, combination.get('t_list')))
                t_list.sort()
                p_list = list(map(int, combination.get('p_list')))
                p_list.sort()
                n_mol_ratio = combination.get('n_mol_ratio')
                if len(smiles_list) > len(set(smiles_list)):
                    raise Exception('Duplicated smiles %s' % smiles_list)

                if not (len(smiles_list) == len(name_list) == len(n_mol_ratio) ):
                    raise Exception('Length of smiles names and n_mol_ratio must be the same')  # Provide name for no or all molecules
                elif not re.match('^[A-Za-z0-9-]+$', ''.join(name_list)):
                    raise Exception('Only letters, numbers and - are valid for names')

                if procedure not in current_app.config['ALLOWED_PROCEDURES']:
                    raise Exception('Invalid procedure for this config')

                smiles_list = [get_canonical_smiles(smiles) for smiles in smiles_list]
                tasks = Task.query.filter(Task.smiles_list == json.dumps(smiles_list)).filter(
                    Task.procedure == procedure). \
                    filter(Task.n_mol_ratio == json.dumps(n_mol_ratio))
                # get t p list from prior task
                if prior_procedure is not None:
                    prior_tasks = Task.query.filter(Task.smiles_list == json.dumps(smiles_list)).filter(
                        Task.procedure == prior_procedure). \
                        filter(Task.n_mol_ratio == json.dumps(n_mol_ratio))
                    if prior_tasks.count() == 1:
                        prior_task = prior_tasks.first()
                    else:
                        continue
                    if t_list == p_list == []:
                        t_list = json.loads(prior_task.t_list)
                        p_list = json.loads(prior_task.p_list)
                    else:
                        for t in t_list:
                            if t not in json.loads(prior_task.t_list):
                                t_list.remove(t)
                        for p in p_list:
                            if p not in json.loads(prior_task.p_list):
                                p_list.remove(p)
                # alter exist task t p list
                if tasks.count() == 1 and alter_task:
                    task = tasks.first()
                    t_list_new = list(set(json.loads(task.t_list) + t_list))
                    t_list_new.sort()
                    p_list_new = list(set(json.loads(task.p_list) + p_list))
                    p_list_new.sort()
                    if t_list_new != json.loads(task.t_list) or p_list_new != json.loads(task.p_list):
                        task.t_list = json.dumps(t_list_new)
                        task.p_list = json.dumps(p_list_new)
                        M += 1
                        try:
                            task.insert_jobs()
                            for job in task.jobs.filter(Job.status == Compute.Status.STARTED):
                                commands = job.prepare()
                                task.commands = json.dumps(commands)
                        except Exception as e:
                            current_app.logger.error('Alter task failed %s: %s' % (self, repr(e)))
                            traceback.print_exc()
                            task.status = Compute.Status.FAILED
                        else:
                            task.stage = Compute.Stage.BUILDING
                            task.status = Compute.Status.DONE
                        db.session.commit()
                # Ignore existed task
                # add a new task
                elif tasks.count() == 0 and create_task:
                    task = Task()
                    if Task.query.count() != 0:
                        task.id = Task.query.order_by(Task.id)[-1].id + 1
                    task.compute_id = self.id
                    task.n_components = len(smiles_list)
                    task.smiles_list = json.dumps(smiles_list)
                    task.procedure = procedure
                    task.t_list = json.dumps(t_list)
                    task.p_list = json.dumps(p_list)
                    task.n_mol_ratio = json.dumps(n_mol_ratio)
                    if name_list is not None:
                        task.name = '_'.join(name_list) + '_' + random_string(4)
                    db.session.add(task)
                    db.session.commit()
                    N += 1

            db.session.commit()
            current_app.logger.info('%i tasks created, %i tasks altered, %s' % (N, M, self))
        except Exception as e:
            db.session.rollback()
            current_app.logger.error('Create tasks failed %s %s' % (self, repr(e)))
            raise Exception('Create tasks failed: ' + repr(e))

    class Stage:
        SUBMITTED = 0
        BUILDING = 1
        RUNNING = 2

        text = {
            SUBMITTED: 'Submitted',
            BUILDING : 'Building',
            RUNNING  : 'Running',
        }

    class Status:
        STARTED = 1
        DONE = 9
        FAILED = -1
        ANALYZED = 10

        text = {
            STARTED : 'Started',
            DONE    : 'Done',
            FAILED  : 'Failed',
            ANALYZED: 'Analyzed',
        }


class Task(db.Model):

    '''
    A Task record corresponds to one molecules and a series of physical states
    A Task record can produce a lots of Job records
    '''
    __tablename__ = 'task'
    id = NotNullColumn(Integer, primary_key=True)
    compute_id = NotNullColumn(Integer, ForeignKey(Compute.id))
    n_components = NotNullColumn(Integer)
    smiles_list = NotNullColumn(Text)
    n_mol_list = Column(Text, nullable=True)
    procedure = NotNullColumn(String(200))
    t_list = Column(Text)
    p_list = Column(Text)
    n_mol_ratio = Column(Text)
    name = NotNullColumn(String(200), default=random_string)
    stage = NotNullColumn(Integer, default=Compute.Stage.SUBMITTED)
    status = NotNullColumn(Integer, default=Compute.Status.DONE)
    commands = Column(Text, nullable=True)
    remark = Column(Text, nullable=True)
    post_result = Column(Text, nullable=True)
    atom_type = Column(Text)

    compute = db.relationship(Compute)
    jobs = db.relationship('Job', lazy='dynamic')

    def __repr__(self):
        return '<Task: %i: %s %s %s>' % (self.id, self.procedure, self.name, self.smiles_list)

    @property
    def dir(self) -> str:
        return os.path.join(Config.WORK_DIR, self.procedure, self.name)

    @property
    def remote_dir(self) -> str:
        if current_app.jobmanager.is_remote:
            return os.path.join(current_app.jobmanager.remote_dir, self.procedure, self.name)
        else:
            return self.dir

    @property
    def prior_task(self):
        def is_subset(a, b):# if a is a subset of b, return True
            for x in a:
                if x not in b:
                    return False
            return True
        prior_procedure = Procedure.prior.get(self.procedure)
        if prior_procedure is None:
            return None

        # TODO Herein we suppose there is no duplicated task. Do not consider T and P
        tasks = Task.query.filter(Task.smiles_list == self.smiles_list).filter(Task.procedure == prior_procedure)
        for task in tasks:
            if is_subset(json.loads(self.t_list), json.loads(task.t_list)) and is_subset(json.loads(self.p_list), json.loads(task.p_list)):
                return task
        return None

    @property
    def n_mol_total(self) -> int:
        return sum(json.loads(self.n_mol_list))

    def get_charge_list(self):
        smiles_list = self.get_smiles_list()
        charge_list = []
        if smiles_list == None:
            return None
        else:
            import pybel
            for smiles in smiles_list:
                mol = pybel.readstring('smi', smiles)
                charge_list.append(mol.charge * Config.CHARGE_SCALE)
        return charge_list

    def get_smiles_list(self):
        if self.smiles_list is None:
            return None
        else:
            return json.loads(self.smiles_list)

    def get_distinct_mol_name_list(self):
        name_list = []
        n = len(json.loads(self.n_mol_list))
        for i in range(n):
            name_list.append('MO%i' % (i))
        return name_list

    def get_LJ_atom_type(self):
        if self.prior_task is None:
            ppf = os.path.join(self.dir, 'build', 'ff.ppf')
        else:
            ppf = os.path.join(self.prior_task.dir, 'build', 'ff.ppf')
        if os.path.exists(ppf):
            LJ_atom_type = []
            for line in open(ppf, 'r').readlines():
                if line.startswith('N12_6:'):
                    atom_type = line.split()[1].split(':')[0]
                    if atom_type not in LJ_atom_type:
                        LJ_atom_type.append(atom_type)
                    else:
                        return
            LJ_atom_type.sort()
            self.atom_type = json.dumps(LJ_atom_type)
            db.session.commit()
        else:
            return

    def get_mol_list(self):
        if self.smiles_list is None:
            return None
        else:
            import pybel
            return [pybel.readstring('smi', smiles) for smiles in json.loads(self.smiles_list)]

    def get_post_result(self):
        if self.post_result is None:
            return None
        else:
            return json.loads(self.post_result)

    def get_t_list(self):
        return json.loads(self.t_list)

    def get_p_list(self):
        if self.p_list != '[]':
            return json.loads(self.p_list)
        else:
            return [None]

    @classmethod
    def get_task(cls, smiles, procedure='npt'):
        return Task.query.filter(Task.smiles_list == json.dumps([smiles])).filter(Task.procedure == procedure).first()

    def build(self):
        current_app.logger.info('Build %s' % self)
        try:
            cd_or_create_and_cd(self.dir)
            cd_or_create_and_cd('build')

            self.stage = Compute.Stage.BUILDING
            self.status = Compute.Status.STARTED
            db.session.commit()

            sim = init_simulation(self.procedure)
            smiles_list = json.loads(self.smiles_list)
            if self.n_mol_list is None:
                if self.procedure in ['npt', 'npt-multi', 'npt-v-rescale']:
                    sim.set_system(smiles_list, n_atoms=Config.NATOMS, n_mols=Config.NMOLS)
                    sim.build(ppf=Config.PPF)
                    self.n_mol_list = json.dumps(sim.n_mol_list)
                elif self.procedure == 'nvt-slab':
                    sim.set_system(smiles_list)
                    sim.build(ppf=Config.PPF)
                    self.n_mol_list = json.dumps(sim.n_mol_list)
                elif self.procedure == 'ppm':
                    n_mol_list = json.loads(self.prior_task.n_mol_list)
                    n_mol_list = (np.array(n_mol_list) * 2).tolist()
                    self.n_mol_list = json.dumps(n_mol_list)
                elif self.procedure in ['nvt-multi', 'nvt-multi-2', 'nvt-multi-3']:
                    self.n_mol_list = self.prior_task.n_mol_list
                elif self.procedure == 'npt-2':
                    n_mol_list = json.loads(self.prior_task.n_mol_list)
                    n_mol_list = (np.array(n_mol_list) * 2).tolist()
                    self.n_mol_list = json.dumps(n_mol_list)
                    sim.set_system(smiles_list, n_mol_list=json.loads(self.n_mol_list))
                    sim.build(ppf=Config.PPF)
                elif self.procedure == 'npt-3':
                    n_mol_list = json.loads(self.prior_task.n_mol_list)
                    n_mol_list = (np.array(n_mol_list) * 3).tolist()
                    self.n_mol_list = json.dumps(n_mol_list)
                    sim.set_system(smiles_list, n_mol_list=json.loads(self.n_mol_list))
                    sim.build(ppf=Config.PPF)
            else:
                if self.procedure == 'npt':
                    sim.set_system(smiles_list, n_mol_list=json.loads(self.n_mol_list))
                    sim.build(ppf=Config.PPF)
                else:
                    raise Exception('for procedure=%s, n_mol_list cannot assigned' % (self.procedure))
        except Exception as e:
            current_app.logger.error('Build task failed %s: %s' % (self, repr(e)))
            traceback.print_exc()
            self.status = Compute.Status.FAILED
            db.session.commit()
        else:
            try:
                self.insert_jobs()
                for job in self.jobs.filter(Job.status == Compute.Status.STARTED):
                    commands = job.prepare()
                    self.commands = json.dumps(commands)
            except Exception as e:
                current_app.logger.error('Build task failed %s: %s' % (self, repr(e)))
                traceback.print_exc()
                self.status = Compute.Status.FAILED
            else:
                self.status = Compute.Status.DONE
            db.session.commit()

    def _rebuild(self):
        '''
        For debug
        '''
        cd_or_create_and_cd(self.dir)
        cd_or_create_and_cd('build')

        sim = init_simulation(self.procedure)
        sim.set_system(json.loads(self.smiles_list), n_mol_list=json.loads(self.n_mol_list))
        sim.build()

    def insert_jobs(self):
        if not os.path.exists(os.path.join(self.dir, 'build')):
            raise Exception('Should build simulation box first')

        if json.loads(self.p_list) == []:
            p_list = [None]
        else:
            p_list = json.loads(self.p_list)
        if Job.query.count() != 0:
            job_id = Job.query.order_by(Job.id)[-1].id + 1
        else:
            job_id = 1
        for p in p_list:
            for t in json.loads(self.t_list):
                for i in range(current_app.config['REPEAT_NUMBER']):
                    if Job.query.filter(Job.task_id == self.id).filter(Job.t == t).filter(
                            Job.p == p).filter(Job.repeat_id == i + 1).first() is not None:
                        continue
                    job = Job()
                    job.id = job_id
                    job_id += 1
                    job.task_id = self.id
                    job.t = t
                    job.p = p
                    job.repeat_id = i + 1
                    job.name = '%s-%i-%i' % (self.name, t or 0, p or 0)
                    prior_procedure = Procedure.prior.get(self.procedure)
                    if prior_procedure != None:
                        if self.prior_task == None or self.prior_task.status != Compute.Status.ANALYZED:
                            self.reset()
                            return
                            # job.status = Compute.Status.FAILED
                            # job.result = {'failed': True,
                            # 'reason': 'prior %s job need to be done first' % (prior_procedure)}
                        elif self.prior_task.jobs.filter(and_(Job.t == t, Job.p == p)).first() == None:
                            job.status = Compute.Status.FAILED
                            job.result = {'failed': True,
                                          'reason': 'prior %s job need to be done first in condition t=%f p=%f' % (
                                              prior_procedure, t, p)}
                        elif not self.prior_task.jobs.filter(and_(Job.t == t, Job.p == p)).first().converged:
                            job.status = Compute.Status.FAILED
                            job.result = {'failed': True,
                                          'reason': 'prior %s job do not converge in condition t=%f p=%f' % (
                                              prior_procedure, t, p)}
                    db.session.add(job)
                    db.session.commit()


        db.session.commit()

    # this function is used to extend a task with more repeated jobs using different initial random number
    def check_task_jobs(self):
        current_app.logger.info('check %s' % self)
        if self.jobs.count() != (current_app.config['REPEAT_NUMBER'] * len(self.get_t_list()) * len(self.get_p_list())) and self.jobs.count() != 0:
            current_app.logger.info('add repeated jobs %s' % self)
            self.insert_jobs()
            for job in self.jobs:
                if not os.path.exists(job.dir):
                    job.prepare()
            self.stage = Compute.Stage.BUILDING
            self.status = Compute.Status.DONE
            db.session.commit()

    def get_exist_pbs_number(self):
        i = 0
        while os.path.exists(os.path.join(self.dir, '_job.run-%i.sh' % (i))):
            i += 1
        return i

    def run(self, ignore_pbs_limit=False) -> int:
        '''
        :return: -1  if jobmanager not working or PBS_NJOB_LIMIT reached
                  0  if no job submit or failed
                  N  if N pbs jobs submit
        '''
        if not (self.stage == Compute.Stage.BUILDING and self.status == Compute.Status.DONE):
            raise Exception('Incorrect stage or status: %s %s' %
                            (Compute.Stage.text[self.stage], Compute.Status.text[self.status]))

        current_app.logger.info('Run %s' % self)
        if not current_app.jobmanager.is_working():
            current_app.logger.warning('JobManager not working')
            return -1

        jobs_to_run = self.jobs.filter(Job.status == Compute.Status.STARTED).filter(Job.cycle == 0).filter(Job.pbs_jobs_id == None)
        n_pbs_run = jobs_to_run.count()
        if current_app.config['GMX_MULTI']:
            multi_njob = current_app.config['GMX_MULTI_NJOB']
            n_pbs_run = int(math.ceil(n_pbs_run / multi_njob))

        if not ignore_pbs_limit and \
                current_app.jobmanager.n_running_jobs + n_pbs_run > current_app.config['PBS_NJOB_LIMIT']:
            current_app.logger.warning('PBS_NJOB_LIMIT reached')
            return -1

        os.chdir(self.dir)

        # TODO Be careful, this is confusing. Upload task to remote server. May failed
        if current_app.jobmanager.is_remote:
            if not current_app.jobmanager.upload(remote_dir=self.remote_dir):
                current_app.logger.warning('Failed upload %s' % self)
                return 0

        self.stage = Compute.Stage.RUNNING
        try:
            if not current_app.config['GMX_MULTI']:
                for job in jobs_to_run:
                    job.run()
                    time.sleep(0.2)
            else:
                # Generate sh for multi simulation based on self.commands
                multi_dirs = [job.dir for job in jobs_to_run]
                multi_cmds = json.loads(self.commands)
                
                commands_list = GMX.generate_gpu_multidir_cmds(multi_dirs, multi_cmds,
                                                               n_parallel=multi_njob,
                                                               n_gpu=current_app.jobmanager.ngpu,
                                                               n_omp=current_app.config['GMX_MULTI_NOMP'])
                njobs_command = []
                n = len(multi_dirs)
                while True:
                    n -= multi_njob
                    if n > 0:
                        njobs_command.append(multi_njob)
                    else:
                        njobs_command.append(n + multi_njob)
                        break

                n = self.get_exist_pbs_number()
                for i, commands in enumerate(commands_list):
                    # instead of run directly, we add a record to pbs_job
                    sh = os.path.join(self.dir, '_job.run-%i.sh' % (i + n))
                    pbs_name = '%s-run-%i' % (self.name, (i + n))

                    if current_app.config['GMX_MULTI_NOMP'] is None:
                        n_tasks = None
                    else:
                        n_tasks = current_app.config['GMX_MULTI_NJOB'] * current_app.config['GMX_MULTI_NOMP'] * \
                                  current_app.jobmanager.nprocs_request / current_app.jobmanager.nprocs

                    if current_app.config['PBS_ARGS'][0] != 'gtx':
                        ngpu = None
                    elif njobs_command[i] % 2 == 0:
                        ngpu = 2
                    else:
                        ngpu = 1
                    current_app.jobmanager.generate_sh(self.dir, commands, name=pbs_name, sh=sh, n_tasks=n_tasks, ngpu=ngpu)

                    pbs_job = PbsJob()
                    pbs_job.name = pbs_name
                    pbs_job.sh_file = sh
                    db.session.add(pbs_job)
                    db.session.flush()

                    # save pbs_job_id for jobs
                    # updated jobs will be removed from jobs_to_run
                    for job in jobs_to_run[0: current_app.config['GMX_MULTI_NJOB']]:
                        job.pbs_jobs_id = json.dumps([pbs_job.id])
                    db.session.commit()

                    # submit job, record if success or failed
                    pbs_job.submit(remote_dir=self.remote_dir, local_dir=self.dir)
                    time.sleep(0.2)

        except Exception as e:
            current_app.logger.error('Run task failed %s %s' % (self, repr(e)))
            traceback.print_exc()
            self.status = Compute.Status.FAILED
            db.session.commit()
            return 0
        else:
            self.status = Compute.Status.STARTED
            db.session.commit()
            return n_pbs_run

    def check_finished_multiprocessing(self):
        """
        check if all jobs in this tasks are finished
        if finished, analyze the job
        """
        from multiprocessing import Pool

        current_app.logger.info('Check status Multi %s' % self)
        for job in self.jobs.filter(Job.status == Compute.Status.STARTED):
            try:
                job.check_finished()
            except Exception as e:
                current_app.logger.error('Check job status failed %s %s' % (job, repr(e)))
                traceback.print_exc()

        jobs_to_analyze = self.jobs.filter(Job.status == Compute.Status.DONE).all()

        for job in jobs_to_analyze:
            current_app.logger.info('Analyze %s' % job)

        n_process = 8
        n_group = int(math.ceil(len(jobs_to_analyze) / n_process))
        for i in range(n_group):
            job_group = jobs_to_analyze[i * n_process:(i + 1) * n_process]
            with Pool(n_process) as p:
                return_dicts = p.map(wrapper_cls_func,
                                     [(job, 'analyze_multiprocessing', [],
                                       {'job_dir': job.dir, 'job_procedure': job.task.procedure})
                                      for job in job_group])
            for j, job in enumerate(job_group):
                return_dict = return_dicts[j]
                exception = return_dict.pop('exception')
                if exception is not None:
                    current_app.logger.error('Analyze failed %s %s' % (job, exception))
                for k, v in return_dict.items():
                    setattr(job, k, v)

            db.session.commit()

        # Set status as DONE if all jobs are converged or failed
        # Set status as ANALYZED only if all jobs are converged
        _failed = []
        for job in self.jobs:
            # at least one of the job is failed
            if job.status == Compute.Status.FAILED:
                _failed.append(True)
            else:
                _failed.append(False)
            # There is a job is not converged and not failed, need to extend
            if job.converged == False and job.status != Compute.Status.FAILED:
                break
        else:
            # all jobs are failed
            if set(_failed) == {True}:
                self.status = Compute.Status.FAILED
            # all jobs are converged
            elif set(_failed) == {False}:
                self.status = Compute.Status.ANALYZED
            # some jobs failed and some jobs converged
            else:
                self.status = Compute.Status.DONE
            db.session.commit()

    def reset(self, add_mol_list=None):
        for job in self.jobs:
            if job.pbs_jobs_id is not None and job.is_running():
                for pbs_job in job.pbs_jobs():
                    pbs_job.kill()
            db.session.delete(job)
        try:
            shutil.rmtree(self.dir)
        except:
            current_app.logger.warning('Reset task %s Cannot remove folder: %s' % (self, self.dir))

        if add_mol_list is None:
            self.n_mol_list = None
        else:
            n_mol_list = json.loads(self.n_mol_list)
            if len(add_mol_list) != len(n_mol_list):
                raise Exception('add_mol_list must have same length of n_mol_list')
            for i in range(len(add_mol_list)):
                n_mol_list[i] += add_mol_list[i]
            self.n_mol_list = json.dumps(n_mol_list)
        self.stage = Compute.Stage.SUBMITTED
        self.status = Compute.Status.DONE
        self.commands = None
        self.post_result = None
        db.session.commit()

    def remove(self):
        current_app.logger.info('Remove task %s' % self)
        self.reset()
        db.session.delete(self)
        db.session.commit()

    def post_process(self, force=False, t_list=None, p_list=None, repeat_number=None, overwrite=False):
        if self.post_result is not None and not overwrite:
            return
        if not force:
            if not (self.stage == Compute.Stage.RUNNING and self.status in [Compute.Status.DONE, Compute.Status.ANALYZED]):
                return
        T_list = []
        P_list = []
        result_list = []

        if t_list is None:
            t_list = self.get_t_list()
        if p_list is None:
            p_list = self.get_p_list()

        for t in t_list:
            for p in p_list:
                current_app.logger.info('Post-process, t=%i,p=%i task=%s' % (t, p or 0, self))
                jobs = self.jobs.filter(Job.t == t).filter(Job.p == p).filter(Job.converged)
                if jobs.count() == 0:
                    continue
                if repeat_number is not None and repeat_number < jobs.count():
                    jobs = jobs.limit(repeat_number)

                if self.procedure in ['nvt-multi', 'nvt-multi-2', 'nvt-multi-3', 'npt-multi']:

                    viscosity, score1 = self.post_process_GK(t, p, jobs, property='viscosity')
                    electrical_conductivity, score2 = self.post_process_GK(t, p, jobs,
                                                                           property='electrical conductivity')

                    info_dict = {
                        'diffusion constant': {}, # {name: [diff, stderr]}
                        'Nernst-Einstein electrical conductivity': [],
                    }
                    for name in (json.loads(jobs[0].result)).get('diffusion constant').keys():
                        info_dict.get('diffusion constant').update({name: []})
                    for job in jobs:
                        result = json.loads(job.result)
                        info_dict.get('Nernst-Einstein electrical conductivity').append(
                            result.get('Nernst-Einstein electrical conductivity')[0])
                        for name in result.get('diffusion constant').keys():
                            info_dict.get('diffusion constant').get(name).append(
                                result.get('diffusion constant').get(name)[0])
                    econ = np.array(info_dict.get('Nernst-Einstein electrical conductivity'))
                    info_dict['Nernst-Einstein electrical conductivity'] = [econ.mean(), econ.std()]
                    for name in info_dict.get('diffusion constant').keys():
                        diff = np.array(info_dict.get('diffusion constant').get(name))
                        info_dict.get('diffusion constant')[name] = [diff.mean(), diff.std()]

                    # too slow to use green-kubo method to calculate diffusion constants
                    if Config.DIFF_GK:
                        diffusion_constant_score = {}
                        diff, s = self.post_process_GK(t, p, jobs, property='diffusion constant', name='System')
                        diffusion_constant_score.update({'System': [diff, s]})
                        for name in self.get_distinct_mol_name_list():
                            diff, s = self.post_process_GK(t, p, jobs, property='diffusion constant', name=name)
                            diffusion_constant_score.update({name: [diff, s]})
                        info_dict.update({
                            'diffusion constant-gk and score': diffusion_constant_score, # {name: [diff, score]} # time-decomposition method
                        })
                        diffusion_constant_list = {}
                        diffusion_constant_list['System'] = []
                        for name in self.get_distinct_mol_name_list():
                            diffusion_constant_list[name] = []
                        for job in jobs:
                            result = json.loads(job.result)
                            diffusion_constant_list['System'].append(result['diffusion constant-gk']['System'])
                            for name in self.get_distinct_mol_name_list():
                                diffusion_constant_list[name].append(result['diffusion constant-gk'][name])
                        diffusion_constant_stderr = {}
                        diffusion_constant_stderr['System'] = [np.mean(diffusion_constant_list['System']), np.std(diffusion_constant_list['System'])]
                        for name in self.get_distinct_mol_name_list():
                            diffusion_constant_stderr[name] = [np.mean(diffusion_constant_list[name]), np.std(diffusion_constant_list[name])]
                        info_dict.update({
                            'diffusion constant-gk and stderr': diffusion_constant_stderr, # {name: [diff, stderr]} #
                        })

                    info_dict.update({
                        'viscosity': viscosity,
                        'vis_score': score1,
                        'electrical conductivity': electrical_conductivity,
                        'econ_score': score2,
                    })
                    T_list.append(t)
                    P_list.append(p)
                    result_list.append(info_dict)
                else:
                    T_list.append(t)
                    P_list.append(p)
                    result_list.append(json.loads(jobs[0].result))

        sim = init_simulation(self.procedure)
        post_result, post_info = sim.post_process(T_list=T_list, P_list=P_list, result_list=result_list,
                                                  n_mol_list=json.loads(self.n_mol_list))

        if post_result is not None:
            self.post_result = json.dumps(post_result)
            db.session.commit()
        else:
            print(post_info)
        return post_info

    def post_process_GK(self, t, p, jobs, property=None, fit=True, name=None, plot_not_converged=False):
        from mstools.analyzer.acf import get_block_average, Property_dict
        from mstools.analyzer.fitting import polyval, polyfit, ExpConstfit, ExpConstval
        from mstools.analyzer.plot import plot

        def get_property_mean_std_list(jobs, property=None, weight=0.00, scale=1., name=None):
            from mstools.analyzer.acf import get_t_property_list
            t_list = []
            data_list = []
            for job in jobs:
                tl, data = get_t_property_list(dir=job.dir, property=property, weight=weight, name=name)
                if tl is None or data is None:
                    continue
                t_list.append(tl)
                data_list.append(data)
                if not (tl.size == data.size == t_list[0].size == data_list[0].size):
                    t_list.remove(tl)
                    data_list.remove(data)
            # calculate averaged viscosity and its std over many independent parallel simulations
            mean = np.linspace(0, 0, data_list[0].size)  # averaged value of target property
            stderr = np.linspace(0, 0, data_list[0].size)  # stderr
            for data in data_list:
                mean += data
            mean /= len(data_list)
            for data in data_list:
                stderr += (data - mean) ** 2
            stderr = np.sqrt(stderr / (len(data_list) - 1))
            if property == 'electrical conductivity':
                scale = Config.CHARGE_SCALE ** 2
            return t_list[0], mean*scale, stderr*scale

        t_list, mean, stderr = get_property_mean_std_list(jobs, property=property, name=name)

        if fit:
            # fit the std of data using function y(t)=At^b
            _log_t = np.log(t_list)
            _log_stderr = np.log(stderr)
            c1, s1 = polyfit(_log_t, _log_stderr, 1)

            # fit the data using exponential function
            n_block = len([t for t in t_list if t < 1])
            t_list_block = get_block_average(t_list, n_block=n_block)
            if property == 'viscosity':
                bounds = ([-np.inf, 0, 0], [0, np.inf, np.inf])
                c2, s2 = ExpConstfit(t_list_block[2:], get_block_average(mean, n_block=n_block)[2:],
                                 bounds=bounds)
            elif property == 'electrical conductivity':
                # bounds = ([0, 0, 0], [np.inf, np.inf, np.inf])
                bounds = ([0, 0, 0], [100, 100, 100])
                c2, s2 = ExpConstfit(t_list_block[2:], get_block_average(mean, n_block=n_block)[2:],
                                 bounds=bounds)
            elif property == 'diffusion constant':
                factor = math.floor(math.log10(mean.mean()))
                bounds = ([0, 0, 0], [100, 100, 100])
                c2, s2 = ExpConstfit(t_list_block[2:], get_block_average(mean * 10 ** (-factor), n_block=n_block)[2:],
                                 bounds=bounds)
                c2[0] *= 10 ** (factor)
                c2[1] *= 10 ** (factor)
            else:
                return None, 0
            if abs(ExpConstval(t_list[-1], c2) - ExpConstval(t_list[-1] / 2, c2)) / ExpConstval(t_list[-1], c2) > 0.01:
                if plot_not_converged:
                    plot(t_list, mean, ExpConstval(t_list, c2))
                return None, 0.

            t_p_dir = os.path.join(self.dir, '%i-%i' % (t, p))
            if name is not None:
                file_name1 = Property_dict.get(property).get('abbr') + '-%s.txt' % (name)
            else:
                file_name1 = Property_dict.get(property).get('abbr') + '.txt'
            f1 = open(os.path.join(t_p_dir, file_name1), 'w')
            f1.write('#time(ps)\t%s\tfit\n' % (Property_dict.get(property).get('property_unit')))

            if name is not None:
                file_name2 = Property_dict.get(property).get('abbr') + '-%s-stderr.txt' % (name)
            else:
                file_name2 = Property_dict.get(property).get('abbr') + '-stderr.txt'
            f2 = open(os.path.join(t_p_dir, file_name2), 'w')
            f2.write('#time(ps)\t%s_stderr\tfit\n' % (Property_dict.get(property).get('property_unit')))
            for i in range(mean.size):
                f1.write('%#.5e\t%#.5e\t%#.5e\n' % (t_list[i], mean[i], ExpConstval(t_list[i], c2)))
                f2.write('%#.5e\t%#.5e\t%#.5e\n' % (t_list[i], stderr[i], np.exp(polyval(np.log(t_list[i]), c1))))
            return c2[1], s2
        else:
            t_p_dir = os.path.join(self.dir, '%i-%i' % (t, p))
            f1 = open(os.path.join(t_p_dir, '%s.txt' % (self.Property_dict.get(property).get('abbr'))), 'w')
            f1.write('#time(ps)\t%s\tfit\n' % (self.Property_dict.get(property).get('property_unit')))
            f2 = open(os.path.join(t_p_dir, '%s_stderr.txt' % (self.Property_dict.get(property).get('abbr'))), 'w')
            f2.write('#time(ps)\t%s_stderr\tfit\n' % (self.Property_dict.get(property).get('property_unit')))
            for i in range(mean.size):
                f1.write('%f\t%f\n' % (t_list[i], mean[i]))
                f2.write('%f\t%f\n' % (t_list[i], stderr[i]))
            return 0, 0

    def get_post_data(self, T=298, P=1):
        sim = init_simulation(self.procedure)
        post_data = sim.get_post_data(post_result=self.get_post_result(), T=T, P=P, smiles_list=self.get_smiles_list())
        return post_data


class Job(db.Model):
    '''
    A Job record corresponds to one molecule and one physical state (T and P)
    '''
    __tablename__ = 'job'
    id = NotNullColumn(Integer, primary_key=True)
    task_id = NotNullColumn(Integer, ForeignKey(Task.id))
    t = Column(Integer, nullable=True)
    p = Column(Integer, nullable=True)
    time = NotNullColumn(DateTime, default=datetime.now)
    name = NotNullColumn(String(200), default=random_string)
    cycle = NotNullColumn(Integer, default=0)
    status = NotNullColumn(Integer, default=Compute.Status.STARTED)
    converged = NotNullColumn(Boolean, default=False)
    result = Column(Text, nullable=True)
    pbs_jobs_id = Column(Text)
    repeat_id = Column(Integer, nullable=True)
    bugfix = NotNullColumn(Boolean, default=False)

    task = db.relationship(Task)

    def __repr__(self):
        return '<Job: %i: %s %i>' % (self.id, self.name, self.cycle)

    @property
    def dir(self) -> str:
        dir_name = '%i-%i' % (self.t, self.p or 0)
        repeat_name = 'repeat-%i' % (self.repeat_id)
        return os.path.join(self.task.dir, dir_name, repeat_name)

    @property
    def remote_dir(self) -> str:
        dir_name = '%i-%i' % (self.t, self.p or 0)
        repeat_name = 'repeat-%i' % (self.repeat_id)
        return os.path.join(self.task.remote_dir, dir_name, repeat_name)

    @property
    def prior_job(self):
        prior_task = self.task.prior_task
        if prior_task is None:
            return None

        return prior_task.jobs.filter(
                and_(Job.t == self.t,
                     Job.p == self.p
                     )
        ).first()

    @property
    def need_extend(self) -> bool:
        return self.status == Compute.Status.ANALYZED and not self.converged

    @property
    def get_nmol_list(self):
        return json.loads(self.task.n_mol_list)

    def get_charge_list(self):
        return self.task.get_charge_list()

    def is_charged_system(self):
        if set(self.get_charge_list()) == {0}:
            return False
        else:
            return True

    def get_result(self):
        if self.result is None:
            return None
        else:
            return json.loads(self.result)

    def get_simulation_part(self, name): # The job with same value of this function can be simulated together using GPU.
        if name == 'nvt-slab':
            name = 'nvt'
        log = os.path.join(self.dir, '%s.log' % name)
        f = open(log, 'r')
        n = 1
        for line in f.readlines():
            if line.startswith('Started mdrun'):
                n += 1
        return n

    def get_vis_with_weight(self, weight=None):
        if weight is None:
            return
        f_acf = os.path.join(self.dir, 'P_acf.txt')
        if not os.path.exists(f_acf):
            return
        acf_info = pd.read_csv(f_acf, sep='\s+', header=0)
        t_list = np.array(acf_info['#time(ps)'])
        acf_list = np.array(acf_info['ACF(Pab)'])
        dt = t_list[1] - t_list[0]
        temperature = self.t
        gro = os.path.join(self.dir, 'nvt.gro')
        f = open(gro, 'r')
        box = f.readlines()[-1].split()
        volume = float(box[0]) * float(box[1]) * float(box[2])
        convert = 6.022 * 0.001 * volume / (8.314 * temperature)
        vis = convert * acf_list[0] * dt / 2
        f_vis = os.path.join(self.dir, 'vis-%.2f.txt' % (weight))
        f = open(f_vis, 'w')
        f.write('#time(ps)\tviscosity(mPa·s)\n')
        for i in range(1, len(t_list)):
            f.write('%f\t%f\n' % (t_list[i]-0.5*dt, vis))
            if t_list[i] <= 1:
                vis += convert * acf_list[i] * dt
            else:
                vis += convert * acf_list[i] * dt * t_list[i] ** (-weight)

    def prepare(self):
        if self.status == Compute.Status.STARTED:
            prior_job = self.prior_job
            if prior_job is None:
                prior_job_dir = None
            else:
                if not prior_job.converged:
                    current_app.logger.warning('Prepare job %s Prior job not converged' % repr(self))
                    prior_job_dir = None
                else:
                    prior_job_dir = prior_job.dir

            cd_or_create_and_cd(self.dir)
            sim = init_simulation(self.task.procedure)

            # Using Temperature dependent parameters or not
            if Config.DFF_TABLE == 'IL':
                T_basic = 350
                drde = True
            elif Config.DFF_TABLE == 'MGI':
                T_basic = 298
                drde = True
            else:
                return
            # prepare
            if self.task.procedure in ['npt', 'npt-2', 'npt-3']:
                commands = sim.prepare(model_dir='../../build', T=self.t, P=self.p, jobname=self.name,
                                       drde=drde, T_basic=T_basic)
            elif self.task.procedure == 'npt-v-rescale':
                commands = sim.prepare(model_dir='../../build', T=self.t, P=self.p, jobname=self.name,
                                       dt=0.001, nst_trr=50, tcoupl='v-rescale', drde=drde,
                                       acf=True, mstools_dir=Config.MS_TOOLS_DIR)
            elif self.task.procedure == 'npt-multi':
                commands = sim.prepare(model_dir='../../build', T=self.t, P=self.p, jobname=self.name, dt=0.001,
                                       nst_trr=50, nst_edr=5, tcoupl='v-rescale', drde=drde, random_seed=self.repeat_id,
                                       acf=True, diff_gk=Config.DIFF_GK, mstools_dir=Config.MS_TOOLS_DIR)
            elif self.task.procedure in ['nvt-multi', 'nvt-multi-2', 'nvt-multi-3']:
                commands = sim.prepare(prior_job_dir=prior_job_dir, T=self.t, jobname=self.name, gro='npt.gro',
                                       tcoupl='v-rescale', random_seed=self.repeat_id, nst_run=int(2E6),
                                       acf=True, diff_gk=Config.DIFF_GK, mstools_dir=Config.MS_TOOLS_DIR)
            elif self.task.procedure == 'ppm':
                commands = sim.prepare(prior_job_dir=prior_job_dir, T=self.t, P=self.p, jobname=self.name,
                                       gro='npt.gro', replicate=(1, 1, 2), random_seed=self.repeat_id)
            elif self.task.procedure == 'nvt-slab':
                commands = sim.prepare(model_dir='../../build', T=self.t, jobname=self.name,
                                       drde=drde)
            else:
                return
            # Using LJ 9-6 potential
            if Config.LJ96:
                sim.gmx._LJ96 = True
                if self.task.procedure in ['npt', 'npt-2', 'npt-3', 'npt-v-rescale', 'npt-multi']:
                    sim.gmx.modify_lj96(['topol.itp'],
                                        ['topol.top', 'topol-hvap.top'],
                                        ['grompp-em.mdp', 'grompp-anneal.mdp', 'grompp-eq.mdp', 'grompp-npt.mdp',
                                         'grompp-hvap.mdp'],
                                        ['em.xvg', 'anneal.xvg', 'eq.xvg', 'npt.xvg', 'hvap.xvg'])
                elif self.task.procedure in ['nvt-multi', 'nvt-multi-2', 'nvt-multi-3']:
                    sim.gmx.modify_lj96([], [], ['grompp-eq.mdp', 'grompp-nvt.mdp'], ['eq.xvg', 'nvt.xvg'])
                elif self.task.procedure == 'ppm':
                    mdp_list = ['grompp-eq-%.3f.mdp' % (a) for a in sim.amplitudes_steps.keys()] \
                               + ['grompp-ppm-%.3f.mdp' % (a) for a in sim.amplitudes_steps.keys()]
                    xvg_list = ['eq-%.3f.xvg' % (a) for a in sim.amplitudes_steps.keys()] \
                               + ['ppm-%.3f.xvg' % (a) for a in sim.amplitudes_steps.keys()]
                    sim.gmx.modify_lj96([], [], mdp_list, xvg_list)

            # all jobs under one task share same commands, save this for GMX -multidir simulation
            return commands

    def run(self) -> bool:
        try:
            os.chdir(self.dir)
        except:
            raise Exception('Should prepare job first')

        # instead of run directly, we add a record to pbs_job
        pbs_job = PbsJob()
        pbs_job.name = self.name
        pbs_job.sh_file = os.path.join(self.dir, current_app.jobmanager.sh)
        db.session.add(pbs_job)
        db.session.flush()

        self.pbs_jobs_id = json.dumps([pbs_job.id])
        self.status = Compute.Status.STARTED
        db.session.commit()

        # TODO Be careful. Upload files to remote server. May failed
        if current_app.jobmanager.is_remote:
            if not current_app.jobmanager.upload(remote_dir=self.remote_dir):
                current_app.logger.warning('Failed upload %s' % self)
                return False

        # submit job, record if success or failed
        pbs_job.submit(remote_dir=self.remote_dir, local_dir=self.dir)
        return pbs_job.submitted

    def _rerun(self):
        '''
        For debug
        '''

        if self.task.procedure != 'npt':
            raise ('Invalid')

        self.task.status = 1
        self.status = 1
        self.cycle = 0
        self.converged = False
        self.result = None
        db.session.commit()

        cd_or_create_and_cd(self.dir)

        commands = []

        T = self.t
        P = self.p
        dt = 0.002
        top = 'topol.top'
        nst_run = int(5E5)
        nst_edr = 100
        nst_trr = int(5E4)
        nst_xtc = 1000
        jobname = self.name

        sim = init_simulation(self.task.procedure)

        nprocs = sim.jobmanager.nprocs
        # NPT production with Langevin thermostat and Parrinello-Rahman barostat
        sim.gmx.prepare_mdp_from_template('t_npt.mdp', mdp_out='grompp-npt.mdp', T=T, P=P,
                                          dt=dt, nsteps=nst_run, nstenergy=nst_edr, nstxout=nst_trr, nstvout=nst_trr,
                                          nstxtcout=nst_xtc, restart=True)
        cmd = sim.gmx.grompp(mdp='grompp-npt.mdp', gro='npt.gro', top=top, tpr_out='npt.tpr',
                             cpt='npt.cpt', get_cmd=True)
        commands.append(cmd)
        cmd = sim.gmx.mdrun(name='npt', nprocs=nprocs, get_cmd=True)
        commands.append(cmd)

        # Rerun enthalpy of vaporization
        commands.append('export GMX_MAXCONSTRWARN=-1')

        top_hvap = 'topol-hvap.top'
        sim.gmx.prepare_mdp_from_template('t_npt.mdp', mdp_out='grompp-hvap.mdp', nstxtcout=0, restart=True)
        cmd = sim.gmx.grompp(mdp='grompp-hvap.mdp', gro='npt.gro', top=top_hvap, tpr_out='hvap.tpr', get_cmd=True)
        commands.append(cmd)
        # Use OpenMP instead of MPI when rerun hvap
        cmd = sim.gmx.mdrun(name='hvap', nprocs=nprocs, n_omp=nprocs, rerun='npt.xtc', get_cmd=True)
        commands.append(cmd)

        sim.jobmanager.generate_sh(os.getcwd(), commands, name=jobname or self.procedure)

        self.run()

    def extend(self, hipri=False):
        if not self.need_extend:
            current_app.logger.warning('Will not extend %s Status: %s' % (self, Compute.Status.text[self.status]))
            return

        if self.cycle >= Config.EXTEND_CYCLE_LIMIT:
            current_app.logger.warning('Will not extend %s EXTEND_CYCLE_LIMIT reached' % self)
            return

        os.chdir(self.dir)

        pbs_name = '%s-repeat%i-%i' % (self.name, self.repeat_id, self.cycle + 1)
        sh = os.path.join(self.dir, '_job.extend-%i.sh' % (self.cycle + 1))

        # instead of run directly, we add a record to pbs_job

        sim = init_simulation(self.task.procedure, extend=True)
        if Config.LJ96:
            sim.gmx._LJ96 = True
        info = json.loads(self.result)
        if set(info.get('continue')) != [False]:
            sim.extend(jobname=pbs_name, sh=sh, info=info, dt=sim.dt, hipri=hipri)

        pbs_job = PbsJob(extend=True)
        pbs_job.name = pbs_name
        pbs_job.sh_file = sh
        db.session.add(pbs_job)
        db.session.flush()

        self.pbs_jobs_id = json.dumps([pbs_job.id])
        self.cycle += 1
        self.status = Compute.Status.STARTED
        db.session.commit()

        # submit job, record if success or failed
        pbs_job.submit()
        os.chdir(current_app.config['RUN_DIR'])

    def continue_terminated_job(self):
        if self.bugfix:
            return
        if self.result is None:
            return
        if self.task.procedure != 'ppm':
            return
        if current_app.config['BUGFIX_PBS_ARGS'][0] == 'gtx' or current_app.config['BUGFIX_GMX_MULTI']:
            return
        if self.is_running():
            return

        os.chdir(self.dir)
        result = json.loads(self.result)
        name_list = result.get('name')
        length_list = result.get('length') if type(result.get('length')) == list else [result.get('length')]

        for i, name in enumerate(name_list):
            if not os.path.exists('%s.log' % (name)) \
                    or not os.path.exists('%s.cpt' % (name)) \
                    or not get_last_line('%s.log' % (name)).startswith('Finished mdrun')\
                    or 'Fatal error:\n' in open('%s.log' % (name), 'r').readlines():
                break
            if (length_list[i] is None or int(length_list[i]) % 100 != 0) and os.path.exists('%s.log' % (name)) and os.path.exists('%s.cpt' % (name)):
                sim = init_simulation(self.task.procedure, bugfix=True)
                pbs_name = '%s-repeat%i-bugfix' % (self.name, self.repeat_id)
                commands = sim.extend_single(jobname=pbs_name, name=name)
                n_amplitude = len(sim.amplitudes_steps.keys())
                commands += json.loads(self.task.commands)[4*(i+1):4*n_amplitude]
                sh = '_job.bugfix.sh'
                sim.jobmanager.generate_sh(os.getcwd(), commands, name=pbs_name or self.task.procedure, sh=sh)
                pbs_job = PbsJob(bugfix=True)
                pbs_job.name = pbs_name
                pbs_job.sh_file = sh
                db.session.add(pbs_job)
                db.session.flush()
                self.result = None
                self.pbs_jobs_id = json.dumps([pbs_job.id])
                self.status = Compute.Status.STARTED
                self.task.status = Compute.Status.STARTED
                db.session.commit()
                pbs_job.submit()
                break

        self.bugfix = True
        db.session.commit()

    def pbs_jobs(self):
        pbs_jobs_id = json.loads(self.pbs_jobs_id)
        pbs_jobs = PbsJob.query.filter(PbsJob.id.in_(pbs_jobs_id))
        return pbs_jobs

    def is_pbs_generated(self) -> bool:
        if self.pbs_jobs_id is None:
            return False
        else:
            return True

    def is_running(self) -> bool:
        '''
        check if the job is running

        If the pbs job has not been generated, return False
        If the all pbs jobs are not running, return False
        Otherwise, return True
        :return:
        '''
        if not self.is_pbs_generated():
            return False
        pbs_jobs_id = json.loads(self.pbs_jobs_id)
        pbs_jobs = PbsJob.query.filter(PbsJob.id.in_(pbs_jobs_id))
        for pbs_job in pbs_jobs:
            if pbs_job.is_running():
                return True
        return False

    def is_running_finished(self) -> bool:
        if self.is_pbs_generated() and not self.is_running():
            return True
        else:
            return False

    def check_finished(self) -> bool:
        if self.status in (Compute.Status.FAILED, Compute.Status.ANALYZED):
            return True
        elif self.status == Compute.Status.DONE:
            if self.task.Stage == Compute.Stage.RUNNING:
                return True
            else:
                return False
        if not self.is_running_finished():
            return False

        os.chdir(self.dir)

        # TODO Be careful. Download job result from remote server. Not for extending
        # TODO If the download failed, treat it as not finished
        # TODO Currently disabled. Because all files were downloaded manually
        # if self.cycle == 0 and current_app.jobmanager.is_remote:
        #     if not current_app.jobmanager.download(remote_dir=self.remote_dir):
        #         return False

        sim = init_simulation(self.task.procedure)
        if sim.check_finished():
            self.status = Compute.Status.DONE
        else:
            self.status = Compute.Status.FAILED
            current_app.logger.error('Job failed %s' % self)
        db.session.commit()
        return True

    def analyze_multiprocessing(self, job_dir, job_procedure, **kwargs):
        _exception = None
        _status = self.status
        _converged = self.converged
        _result = self.result

        os.chdir(job_dir)
        sim = init_simulation(job_procedure)
        try:
            if self.task.procedure in ['npt-multi', 'nvt-multi', 'nvt-multi-2', 'nvt-multi-3']:
                result = sim.analyze_acf(mstools_dir=Config.MS_TOOLS_DIR, n_mol_list=json.loads(self.task.n_mol_list),
                                         charge_list=self.task.get_charge_list(), current=self.is_charged_system(),
                                         delete_trr=not Config.DEBUG, diff_gk=Config.DIFF_GK)
            else:
                result = sim.analyze(**kwargs)
        except Exception as e:
            # One kine of fail -- simulation fail
            _status = Compute.Status.FAILED
            _exception = repr(e)
        else:
            if self.task.procedure == 'ppm':
                from collections import Counter
                if Counter(result.get('failed'))[True] <= 2 or result.get('failed') == [False, False, False, True, True, True]:
                    if result.get('converged'):
                        # simulation converged and normally ended
                        _status = Compute.Status.ANALYZED
                        _converged = True
                        _result = json.dumps(result)
                        # Clean intermediate files if converged
                        sim.clean()
                    else:
                        # simulation not converged, need to extend simulation
                        _status = Compute.Status.ANALYZED
                        _converged = False
                        _result = json.dumps(result)
                else:
                    # simulation failed, cannot extend
                    _status = Compute.Status.FAILED
                    _result = json.dumps(result)
            else:
                if not result.get('failed')[0]:
                    if not result.get('continue')[0]:
                        # simulation converged and normally ended
                        if self.task.procedure == 'npt-v-rescale':
                            result.update(sim.analyze_diff(n_mol_list=json.loads(self.task.n_mol_list), charge_list=self.task.get_charge_list()))
                        _status = Compute.Status.ANALYZED
                        _converged = True
                        _result = json.dumps(result)
                        # Clean intermediate files if converged
                        sim.clean()
                    else:
                        # simulation not converged, need to extend simulation
                        _status = Compute.Status.ANALYZED
                        _converged = False
                        _result = json.dumps(result)
                else:
                    # simulation failed, cannot extend
                    _status = Compute.Status.FAILED
                    _result = json.dumps(result)


        return {
            'exception': _exception,
            'status': _status,
            'converged': _converged,
            'result': _result,
        }
    '''
    def analyze_simple(self, job_dir, job_procedure, weight=None):
        os.chdir(job_dir)
        sim = init_simulation(job_procedure)
        if self.task.procedure == 'nvt-multi' and weight is not None:
            for w in weight:
                sim.analyze(skip=1, current=self.is_charged_system(), mstools_dir=Config.MS_TOOLS_DIR,temperature=self.t, weight=w)
    '''
    def remove(self):
        if self.pbs_job_id is not None and self.is_running():
            self.pbs_job.kill()

        try:
            shutil.rmtree(self.dir)
        except:
            current_app.logger.warning('Remove job %s Cannot remove folder: %s' % (self, self.dir))

        db.session.delete(self)
        db.session.commit()


class ExtendJob:
    def __init__(self, job, continue_n):
        self.job = job
        self.continue_n = continue_n

    def __le__(self, other):
        return self.continue_n < other.continue_n

    def __eq__(self, other):
        return self.continue_n == other.continue_n

    def __gt__(self, other):
        return self.continue_n > other.continue_n
