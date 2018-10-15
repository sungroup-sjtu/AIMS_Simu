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

NotNullColumn = partial(Column, nullable=False)


def init_simulation(procedure, extend=False):
    from mstools.simulation import gmx as simulationEngine
    kwargs = {'packmol_bin': Config.PACKMOL_BIN,
              'dff_root'   : Config.DFF_ROOT,
              'dff_table'  : Config.DFF_TABLE,
              'gmx_bin'    : current_app.config['GMX_BIN'],
              'gmx_mdrun'  : current_app.config['GMX_MDRUN'],
              'jobmanager' : current_app.jobmanager
              }
    if extend:
        kwargs.update({
            'gmx_bin'   : current_app.config['EXTEND_GMX_BIN'],
            'gmx_mdrun' : current_app.config['EXTEND_GMX_MDRUN'],
            'jobmanager': current_app.jm_extend
        })

    if procedure == 'npt':
        return simulationEngine.Npt(**kwargs)
    elif procedure == 'nvt':
        return simulationEngine.Nvt(**kwargs)
    elif procedure == 'nvt-slab':
        return simulationEngine.NvtSlab(**kwargs)
    else:
        raise Exception('Unknown simulation procedure')


def wrapper_cls_func(cls_func_args_kwargs):
    cls, func, args, kwargs = cls_func_args_kwargs
    func = getattr(cls, func)
    return func(*args, **kwargs)


class PbsJob(db.Model):
    __tablename__ = 'pbs_job'
    id = NotNullColumn(Integer, primary_key=True)
    name = NotNullColumn(String(200))
    sh_file = NotNullColumn(Text)
    extend = NotNullColumn(Boolean, default=False)
    submitted = NotNullColumn(Boolean, default=False)

    def __repr__(self):
        return '<PbsJob: %i: %s %s>' % (self.id, self.name, self.submitted)

    @property
    def jm(self):
        return current_app.jm_extend if self.extend else current_app.jobmanager

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
        return self.jm.is_running(self.name)


class Compute(db.Model):
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

    def create_tasks(self):
        current_app.logger.info('Create tasks from %s' % self)
        try:
            detail = json.loads(self.json_detail)
            procedures = detail['procedures']
            combinations = detail['combinations']
            t_range_default = detail.get('t')
            p_range_default = detail.get('p')

            # check for procedure and prior
            for p in detail['procedures']:
                if p not in Procedure.choices:
                    raise Exception('Invalid simulation procedure: %s' % p)

                    # prior = Procedure.prior.get(p)
                    # if prior is not None and prior not in procedures:
                    #     procedures = [prior] + procedures  # Put prerequisite at first

            N = 0
            for combination in combinations:
                # check for smiles and name
                smiles_list = combination['smiles']
                name_list = combination.get('names')

                if len(smiles_list) > len(set(smiles_list)):
                    raise Exception('Duplicated smiles %s' % smiles_list)

                if name_list is not None:
                    if len(smiles_list) != len(name_list):
                        raise Exception('Length of names not equal to smiles')  # Provide name for no or all molecules
                    elif not re.match('^[A-Za-z0-9-]+$', ''.join(name_list)):
                        raise Exception('Only letters, numbers and - are valid for names')

                t_range = combination.get('t') or t_range_default
                p_range = combination.get('p') or p_range_default

                if t_range is not None:
                    t_min, t_max = t_range
                    if t_min > t_max:
                        raise Exception('Invalid temperature')
                else:
                    t_min, t_max = None, None

                if p_range is not None:
                    p_min, p_max = p_range
                    if p_min > p_max:
                        raise Exception('Invalid pressure')
                else:
                    p_min, p_max = None, None

                for procedure in procedures:
                    if procedure not in current_app.config['ALLOWED_PROCEDURES']:
                        raise Exception('Invalid procedure for this config')

                    # Ignore existed task
                    if Task.query.filter(Task.smiles_list == json.dumps(smiles_list)) \
                            .filter(Task.procedure == procedure).first() is not None:
                        continue

                    task = Task()
                    task.compute_id = self.id
                    task.n_components = len(smiles_list)
                    task.smiles_list = json.dumps(smiles_list)
                    task.procedure = procedure
                    if procedure in Procedure.T_RELEVANT:
                        if t_min is None or t_max is None:
                            raise Exception('Invalid temperature')
                        task.t_min = t_min
                        task.t_max = t_max
                    else:
                        task.t_min = None
                        task.t_max = None
                    if procedure in Procedure.P_RELEVANT:
                        if p_min is None or p_max is None:
                            raise Exception('Invalid pressure')
                        task.p_min = p_min
                        task.p_max = p_max
                    else:
                        task.p_min = None
                        task.p_max = None
                    if name_list is not None:
                        task.name = '_'.join(name_list) + '_' + random_string(4)
                    db.session.add(task)
                    N += 1
            db.session.commit()
            current_app.logger.info('%i tasks created %s' % (N, self))
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
    __tablename__ = 'task'
    id = NotNullColumn(Integer, primary_key=True)
    compute_id = NotNullColumn(Integer, ForeignKey(Compute.id))
    n_components = NotNullColumn(Integer)
    smiles_list = NotNullColumn(Text)
    n_mol_list = Column(Text, nullable=True)
    procedure = NotNullColumn(String(200))
    t_min = Column(Integer, nullable=True)
    t_max = Column(Integer, nullable=True)
    t_interval = Column(Integer, nullable=True)
    p_min = Column(Integer, nullable=True)
    p_max = Column(Integer, nullable=True)
    name = NotNullColumn(String(200), default=random_string)
    cycle = NotNullColumn(Integer, default=0)
    stage = NotNullColumn(Integer, default=Compute.Stage.SUBMITTED)
    status = NotNullColumn(Integer, default=Compute.Status.DONE)
    commands = Column(Text, nullable=True)
    remark = Column(Text, nullable=True)
    post_result = Column(Text, nullable=True)

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
        prior_procedure = Procedure.prior.get(self.procedure)
        if prior_procedure is None:
            return None

        # TODO Herein we suppose there is no duplicated task. Do not consider T and P
        return Task.query.filter(Task.smiles_list == self.smiles_list).filter(Task.procedure == prior_procedure).first()

    @property
    def n_mol_total(self) -> int:
        return sum(json.loads(self.n_mol_list))

    def get_smiles_list(self):
        if self.smiles_list is None:
            return None
        else:
            return json.loads(self.smiles_list)

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
            sim.set_system(json.loads(self.smiles_list))
            sim.build()
            self.n_mol_list = json.dumps(sim.n_mol_list)
        except Exception as e:
            current_app.logger.error('Build task failed %s: %s' % (self, repr(e)))
            traceback.print_exc()
            self.status = Compute.Status.FAILED
            db.session.commit()
        else:
            try:
                self.insert_jobs()
                for job in self.jobs:
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
        cd_or_create_and_cd(self.dir)
        cd_or_create_and_cd('build')

        sim = init_simulation(self.procedure)
        sim.set_system(json.loads(self.smiles_list), n_mol_list=json.loads(self.n_mol_list))
        sim.build()

    def insert_jobs(self):
        if not os.path.exists(os.path.join(self.dir, 'build')):
            raise Exception('Should build simulation box first')

        if self.t_min is None or self.t_max is None:
            T_list = [None]
        else:
            if self.procedure == 'nvt-slab':
                T_list = get_T_list_VLE_from_range(self.t_min, self.t_max)
            else:
                T_list = get_T_list_from_range(self.t_min, self.t_max)

        if self.p_min is None or self.p_max is None:
            P_list = [None]
        else:
            P_list = [1, 50, 100, 250, 500, 750, 1000]

        for p in P_list:
            for t in T_list:
                if Job.query.filter(Job.task_id == self.id).filter(Job.t == t).filter(Job.p == p).first() is not None:
                    continue

                job = Job()
                job.task_id = self.id
                job.t = t
                job.p = p
                job.name = '%s-%i-%i' % (self.name, t or 0, p or 0)
                db.session.add(job)

        db.session.commit()

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

        n_pbs_run = self.jobs.count()
        if current_app.config['GMX_MULTI']:
            n_pbs_run = int(math.ceil(n_pbs_run / current_app.config['GMX_MULTI_NJOB']))

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
                for job in self.jobs:
                    job.run()
                    time.sleep(0.2)
            else:
                # Generate sh for multi simulation based on self.commands
                multi_dirs = [job.dir for job in self.jobs]
                multi_cmds = json.loads(self.commands)

                commands_list = GMX.generate_gpu_multidir_cmds(multi_dirs, multi_cmds,
                                                               n_parallel=current_app.config['GMX_MULTI_NJOB'],
                                                               n_gpu=current_app.jobmanager.ngpu,
                                                               n_omp=current_app.config['GMX_MULTI_NOMP'])

                for i, commands in enumerate(commands_list):
                    # instead of run directly, we add a record to pbs_job
                    sh = os.path.join(self.dir, '_job.run-%i.sh' % i)
                    pbs_name = '%s-run-%i' % (self.name, i)

                    if current_app.config['GMX_MULTI_NOMP'] is None:
                        n_tasks = None
                    else:
                        n_tasks = current_app.config['GMX_MULTI_NJOB'] * current_app.config['GMX_MULTI_NOMP'] * \
                                  current_app.jobmanager.nprocs_request / current_app.jobmanager.nprocs

                    current_app.jobmanager.generate_sh(self.dir, commands, name=pbs_name, sh=sh, n_tasks=n_tasks)

                    pbs_job = PbsJob()
                    pbs_job.name = pbs_name
                    pbs_job.sh_file = sh
                    db.session.add(pbs_job)
                    db.session.flush()

                    # save pbs_job_id for jobs
                    for job in self.jobs[
                               i * current_app.config['GMX_MULTI_NJOB']:(i + 1) * current_app.config['GMX_MULTI_NJOB']]:
                        job.pbs_job_id = pbs_job.id
                        job.status = Compute.Status.STARTED
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

    def extend(self, ignore_pbs_limit=False) -> int:
        '''
        :return: -1  if EXTEND_PBS_NJOB_LIMIT reached
                  0  if no job submit or failed
                  N  if N pbs jobs submit
        '''
        if not (self.stage == Compute.Stage.RUNNING and self.status == Compute.Status.STARTED):
            raise Exception('Incorrect stage or status: %s %s' %
                            (Compute.Stage.text[self.stage], Compute.Status.text[self.status]))

        jobs_extend = []
        # do not extend the task if some job is running
        for job in self.jobs:
            if job.status == Compute.Status.STARTED:
                return 0
            if job.need_extend and job.cycle < Config.EXTEND_CYCLE_LIMIT:
                jobs_extend.append(job)

        # do not extend the task if no job need extend
        if len(jobs_extend) == 0:
            current_app.logger.warning('No job need extend or EXTEND_CYCLE_LIMIT reached %s' % self)
            return 0

        current_app.logger.info('Extend %s' % self)
        if not current_app.jm_extend.is_working():
            current_app.logger.warning('JobManager not working')
            return -1

        n_pbs_extend = len(jobs_extend)
        if current_app.config['EXTEND_GMX_MULTI']:
            n_pbs_extend = math.ceil(n_pbs_extend / current_app.config['EXTEND_GMX_MULTI_NJOB'])

        if not ignore_pbs_limit and \
                current_app.jm_extend.n_running_jobs + n_pbs_extend > current_app.config['EXTEND_PBS_NJOB_LIMIT']:
            current_app.logger.warning('EXTEND_PBS_NJOB_LIMIT reached')
            return -1

        self.cycle += 1
        db.session.commit()

        try:
            if not current_app.config['EXTEND_GMX_MULTI']:
                for job in jobs_extend:
                    job.extend()
                    time.sleep(0.2)
            else:
                multi_dirs = []
                multi_cmds = []
                sim = init_simulation(self.procedure, extend=True)
                for job in jobs_extend:
                    os.chdir(job.dir)
                    multi_dirs.append(job.dir)
                    multi_cmds = sim.extend(jobname='%s-%i' % (job.name, job.cycle + 1),
                                            sh='_job.extend-%i.sh' % (job.cycle + 1))

                commands_list = GMX.generate_gpu_multidir_cmds(multi_dirs, multi_cmds,
                                                               n_parallel=current_app.config['EXTEND_GMX_MULTI_NJOB'],
                                                               n_gpu=current_app.jm_extend.ngpu,
                                                               n_procs=current_app.jm_extend.nprocs)

                os.chdir(self.dir)
                for i, commands in enumerate(commands_list):
                    sh = os.path.join(self.dir, '_job.extend-%i-%i.sh' % (self.cycle, i))
                    pbs_name = '%s-extend-%i-%i' % (self.name, self.cycle, i)

                    current_app.jm_extend.generate_sh(self.dir, commands, name=pbs_name, sh=sh)

                    # instead of run directly, we add a record to pbs_job
                    pbs_job = PbsJob(extend=True)
                    pbs_job.name = pbs_name
                    pbs_job.sh_file = sh
                    db.session.add(pbs_job)
                    db.session.flush()

                    # save pbs_job_id for jobs
                    for job in jobs_extend[i * current_app.config['EXTEND_GMX_MULTI_NJOB']:(i + 1) * current_app.config[
                        'EXTEND_GMX_MULTI_NJOB']]:
                        job.pbs_job_id = pbs_job.id
                        job.cycle += 1
                        job.status = Compute.Status.STARTED
                    db.session.commit()

                    # submit job, record if success or failed
                    pbs_job.submit()
                    time.sleep(0.2)

        except Exception as e:
            current_app.logger.error('Extend task failed %s %s' % (self, repr(e)))
            traceback.print_exc()
            return 0
        else:
            return n_pbs_extend

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
            for i, job in enumerate(job_group):
                return_dict = return_dicts[i]
                exception = return_dict.pop('exception')
                if exception is not None:
                    current_app.logger.error('Analyze failed %s %s' % (job, exception))
                for k, v in return_dict.items():
                    setattr(job, k, v)

        db.session.commit()

        # Set status as DONE if all jobs are converged or failed
        # Set status as ANALYZED only if all jobs are converged
        _failed = False
        for job in self.jobs:
            if job.status == Compute.Status.FAILED:
                _failed = True
            if not (job.converged or job.status == Compute.Status.FAILED):
                break
        else:
            self.status = Compute.Status.DONE if _failed else Compute.Status.ANALYZED
            db.session.commit()

    def reset(self):
        for job in self.jobs:
            if job.pbs_job_id is not None and job.is_running():
                job.pbs_job.kill()
            db.session.delete(job)
        try:
            shutil.rmtree(self.dir)
        except:
            current_app.logger.warning('Reset task %s Cannot remove folder: %s' % (self, self.dir))

        self.n_mol_list = None
        self.cycle = 0
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

    def post_process(self, force=False):
        if not force:
            if not (self.stage == Compute.Stage.RUNNING and self.status == Compute.Status.ANALYZED):
                return

        T_list = []
        P_list = []
        result_list = []
        for job in self.jobs:
            if not job.converged:
                continue
            T_list.append(job.t)
            P_list.append(job.p)
            result_list.append(json.loads(job.result))

        sim = init_simulation(self.procedure)
        post_result, post_info = sim.post_process(T_list=T_list, P_list=P_list, result_list=result_list)

        if post_result is not None:
            self.post_result = json.dumps(post_result)
            db.session.commit()

        return post_info

    def get_post_data(self, T=298, P=1):
        sim = init_simulation(self.procedure)
        post_data = sim.get_post_data(post_result=json.loads(self.post_result), T=T, P=P,
                                      smiles_list=json.loads(self.smiles_list), n_mol_list=json.loads(self.n_mol_list))
        return post_data


class Job(db.Model):
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
    pbs_job_id = Column(Integer, ForeignKey(PbsJob.id), nullable=True)

    task = db.relationship(Task)
    pbs_job = db.relationship(PbsJob)

    def __repr__(self):
        return '<Job: %i: %s %i>' % (self.id, self.name, self.cycle)

    @property
    def dir(self) -> str:
        dir_name = '%i-%i' % (self.t or 0, self.p or 0)
        return os.path.join(self.task.dir, dir_name)

    @property
    def remote_dir(self) -> str:
        dir_name = '%i-%i' % (self.t or 0, self.p or 0)
        return os.path.join(self.task.remote_dir, dir_name)

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

    def is_running(self) -> bool:
        '''
        If the pbs job missing, return False
        If the pbs job has finished, return False
        Otherwise, return True
        :return:
        '''
        if self.pbs_job_id is None:
            return True
        # Sometimes the pbs_job is missing. I don't know why. Treat it as finished
        elif self.pbs_job is None:
            current_app.logger.warning('%s PbsJob missing' % self)
            return False
        elif not self.pbs_job.submitted:
            return True
        else:
            return self.pbs_job.is_running()

    def get_result(self):
        if self.result is None:
            return None
        else:
            return json.loads(self.result)

    def prepare(self):
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
        commands = sim.prepare(model_dir='../build', T=self.t, P=self.p, jobname=self.name,
                               prior_job_dir=prior_job_dir, drde=True)  # Temperature dependent parameters

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

        self.pbs_job_id = pbs_job.id
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

    def extend(self):
        if not self.need_extend:
            current_app.logger.warning('Will not extend %s Status: %s' % (self, Compute.Status.text[self.status]))
            return

        if self.cycle >= Config.EXTEND_CYCLE_LIMIT:
            current_app.logger.warning('Will not extend %s EXTEND_CYCLE_LIMIT reached' % self)
            return

        os.chdir(self.dir)

        pbs_name = '%s-%i' % (self.name, self.cycle + 1)
        sh = os.path.join(self.dir, '_job.extend-%i.sh' % (self.cycle + 1))

        # instead of run directly, we add a record to pbs_job

        sim = init_simulation(self.task.procedure, extend=True)
        sim.extend(jobname=pbs_name, sh=sh)

        pbs_job = PbsJob(extend=True)
        pbs_job.name = pbs_name
        pbs_job.sh_file = sh
        db.session.add(pbs_job)
        db.session.flush()

        self.pbs_job_id = pbs_job.id
        self.cycle += 1
        self.status = Compute.Status.STARTED
        db.session.commit()

        # submit job, record if success or failed
        pbs_job.submit()

    def check_finished(self) -> bool:
        if self.status in (Compute.Status.DONE, Compute.Status.FAILED, Compute.Status.ANALYZED):
            return True
        if self.is_running():
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

    def analyze(self, rerun=False, **kwargs):
        if self.status == Compute.Status.ANALYZED and not rerun:
            current_app.logger.warning('Will not analyze %s Already analyzed' % self)
            return
        if self.status == Compute.Status.STARTED:
            current_app.logger.warning('Will not analyze %s Job is still running' % self)
            return
        if self.status == Compute.Status.FAILED:
            current_app.logger.warning('Will not analyze %s Job failed' % self)
            return

        os.chdir(self.dir)

        sim = init_simulation(self.task.procedure)
        try:
            result = sim.analyze(**kwargs)
        except:
            # One kine of fail -- simulation fail
            self.status = Compute.Status.FAILED
            db.session.commit()
            raise

        if result is None:
            self.status = Compute.Status.ANALYZED
            self.converged = False
        elif result.get('failed') == True:
            # Another kine of fail -- result not what we want. Record the reason
            self.status = Compute.Status.FAILED
            self.result = json.dumps(result)
        else:
            self.status = Compute.Status.ANALYZED
            self.converged = True
            self.result = json.dumps(result)
            # Clean intermediate files if converged
            sim.clean()

        db.session.commit()

    def analyze_multiprocessing(self, job_dir, job_procedure, **kwargs):
        _exception = None
        _status = self.status
        _converged = self.converged
        _result = self.result

        os.chdir(job_dir)
        sim = init_simulation(job_procedure)
        try:
            result = sim.analyze(**kwargs)
        except Exception as e:
            # One kine of fail -- simulation fail
            _status = Compute.Status.FAILED
            _exception = repr(e)
        else:
            if result is None:
                _status = Compute.Status.ANALYZED
                _converged = False
            elif result.get('failed') == True:
                # Another kine of fail -- result not what we want. Record the reason
                _status = Compute.Status.FAILED
                _result = json.dumps(result)
            else:
                _status = Compute.Status.ANALYZED
                _converged = True
                _result = json.dumps(result)
                # Clean intermediate files if converged
                sim.clean()

        return {'exception': _exception,
                'status'   : _status,
                'converged': _converged,
                'result'   : _result,
                }

    def remove(self):
        if self.pbs_job_id is not None and self.is_running():
            self.pbs_job.kill()

        try:
            shutil.rmtree(self.dir)
        except:
            current_app.logger.warning('Remove job %s Cannot remove folder: %s' % (self, self.dir))

        db.session.delete(self)
        db.session.commit()
