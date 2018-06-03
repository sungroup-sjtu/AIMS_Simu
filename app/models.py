import json
import re
import shutil
import sys
import time
from datetime import datetime
from functools import partial

from sqlalchemy import Column, ForeignKey, Integer, Text, String, Boolean, DateTime, and_

from config import Config
from . import db, log

sys.path.append(Config.MS_TOOLS_DIR)
from mstools.simulation.procedure import Procedure
from mstools.utils import *

NotNullColumn = partial(Column, nullable=False)

from mstools.jobmanager import *
from mstools.wrapper import GMX

if Config.PBS_MANAGER == 'slurm':
    PBS = Slurm
elif Config.PBS_MANAGER == 'torque':
    PBS = Torque
elif Config.PBS_MANAGER == 'local':
    PBS = Local
else:
    raise Exception('Job manager not supported')
jobmanager = PBS(queue_list=Config.PBS_QUEUE_LIST, env_cmd=Config.PBS_ENV_CMD)
jm_extend = PBS(queue_list=Config.EXTEND_PBS_QUEUE_LIST, env_cmd=Config.PBS_ENV_CMD)
if hasattr(Config, 'PBS_SUBMIT_CMD'):
    jobmanager.submit_cmd = Config.PBS_SUBMIT_CMD
    jm_extend.submit_cmd = Config.PBS_SUBMIT_CMD


def init_simulation(procedure, extend=False):
    from mstools.simulation import gmx as simulationEngine
    kwargs = {'packmol_bin': Config.PACKMOL_BIN,
              'dff_root'   : Config.DFF_ROOT,
              'dff_table'  : Config.DFF_TABLE,
              'gmx_bin'    : Config.GMX_BIN,
              'jobmanager' : jobmanager
              }
    if extend:
        kwargs.update({
            'gmx_bin'   : Config.EXTEND_GMX_BIN,
            'jobmanager': jm_extend
        })

    if procedure == 'npt':
        return simulationEngine.Npt(**kwargs)
    elif procedure == 'nvt':
        return simulationEngine.Nvt(**kwargs)
    elif procedure == 'nvt-slab':
        return simulationEngine.NvtSlab(**kwargs)
    else:
        raise Exception('Unknown simulation procedure')


def wrapper_cls_func(cls_funcname_args_kwargs):
    cls, funcname, args, kwargs = cls_funcname_args_kwargs
    func = getattr(cls, funcname)
    return func(*args, **kwargs)


class PbsJob(db.Model):
    __tablename__ = 'pbs_job'
    id = NotNullColumn(Integer, primary_key=True)
    name = NotNullColumn(String(200))
    sh_file = NotNullColumn(Text)
    sh_content = Column(Text, nullable=True)
    submitted = NotNullColumn(Boolean, default=False)

    def __repr__(self):
        return '<PbsJob: %i: %s %s>' % (self.id, self.name, self.submitted)

    def submit(self):
        if jobmanager.submit(self.sh_file):
            self.submitted = True
        else:
            self.submitted = False
            log.warning('Submit PBS job failed %s' % self.sh_file)
        db.session.commit()


class Compute(db.Model):
    __tablename__ = 'compute'
    id = NotNullColumn(Integer, primary_key=True)
    web_id = NotNullColumn(Integer)
    web_user_id = NotNullColumn(Integer)
    web_ip = NotNullColumn(String(200))
    time = NotNullColumn(DateTime, default=datetime.now)
    json = NotNullColumn(Text)
    remark = NotNullColumn(Text)

    tasks = db.relationship('Task', lazy='dynamic')

    def __repr__(self):
        return '<Compute: %i: %s>' % (self.id, self.remark)

    def create_tasks(self):
        log.info('Create tasks from %s' % self)
        try:
            detail = json.loads(self.json)
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
                    # ignore existed task
                    if Task.query.filter(
                            and_(Task.smiles_list == json.dumps(smiles_list),
                                 Task.procedure == procedure
                                 )
                    ).first() is not None:
                        continue

                    task = Task()
                    task.compute_id = self.id
                    task.n_components = len(smiles_list)
                    task.smiles_list = json.dumps(smiles_list)
                    task.procedure = procedure
                    if procedure in Procedure.T_RELEVANT:
                        if t_min == None or t_max == None:
                            raise Exception('Invalid temperature')
                        task.t_min = t_min
                        task.t_max = t_max
                    else:
                        task.t_min = None
                        task.t_max = None
                    if procedure in Procedure.P_RELEVANT:
                        if p_min == None or p_max == None:
                            raise Exception('Invalid pressure')
                        task.p_min = p_min
                        task.p_max = p_max
                    else:
                        task.p_min = None
                        task.p_max = None
                    if name_list is not None:
                        task.name = '_'.join(name_list) + '_' + random_string(4)
                    db.session.add(task)
            db.session.commit()
        except Exception as e:
            db.session.rollback()
            log.error('Create tasks failed %s %s' % (self, repr(e)))
            raise Exception('Create tasks failed: ' + repr(e))

    class Stage:
        SUBMITTED = 0
        BUILDING = 1
        RUNNING = 2

        text = {
            SUBMITTED: 'Submitted',
            BUILDING : 'Building...',
            RUNNING  : 'Running...',
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
        if self.procedure == 'npt':
            return os.path.join(Config.WORK_DIR, self.name)
        else:
            return os.path.join(Config.WORK_DIR, self.procedure, self.name)

    @property
    def prior_task(self):
        prior_procedure = Procedure.prior.get(self.procedure)
        if prior_procedure is None:
            return None

        # TODO Herein we suppose there is no duplicated task. Do not consider T and P
        return Task.query.filter(
                and_(Task.smiles_list == self.smiles_list,
                     Task.procedure == prior_procedure,
                     # Task.t_min == self.t_min,
                     # Task.t_max == self.t_max,
                     # Task.t_interval == self.t_interval,
                     # Task.p_min == self.p_max,
                     # Task.p_max == self.p_max
                     )
        ).first()

    @property
    def n_mol_total(self) -> int:
        n_mol = 0
        for i in json.loads(self.n_mol_list):
            n_mol += i
        return n_mol

    def build(self):
        log.info('Build %s' % self)
        try:
            cd_or_create_and_cd(self.dir)
            cd_or_create_and_cd('build')

            self.stage = Compute.Stage.BUILDING
            self.status = Compute.Status.STARTED
            db.session.commit()

            simulation = init_simulation(self.procedure)
            simulation.set_system(json.loads(self.smiles_list))
            simulation.build()
            self.n_mol_list = json.dumps(simulation.n_mol_list)
        except Exception as e:
            log.error('Build task failed %s: %s' % (self, repr(e)))
            self.status = Compute.Status.FAILED
            db.session.commit()
        else:
            try:
                self.insert_jobs()
                for job in self.jobs:
                    job.prepare()
            except Exception as e:
                log.error('Build task failed %s: %s' % (self, repr(e)))
                self.status = Compute.Status.FAILED
                db.session.commit()
            else:
                self.status = Compute.Status.DONE
                db.session.commit()

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

        try:
            db.session.commit()
        except:
            db.session.rollback()
            raise

    def run(self, ignore_pbs_limit=False) -> bool:
        log.info('Run %s' % self)

        n_pbs_run = self.jobs.count()
        if Config.GMX_MULTI:
            n_pbs_run = math.ceil(n_pbs_run / Config.GMX_MULTI_NJOB)

        if not ignore_pbs_limit and jobmanager.n_running_jobs + n_pbs_run > Config.PBS_NJOB_LIMIT:
            log.warning('PBS_NJOB_LIMIT reached')
            return False

        if self.stage != Compute.Stage.BUILDING:
            print('Incorrect stage: %s' % Compute.Stage.text[self.stage])
            return False
        elif self.status != Compute.Status.DONE:
            print('Incorrect status: %s' % Compute.Status.text[self.status])
            return False

        self.stage = Compute.Stage.RUNNING
        try:
            os.chdir(self.dir)

            if not Config.GMX_MULTI:
                for job in self.jobs:
                    job.run()
                    time.sleep(0.2)
            else:
                multi_dirs = [job.dir for job in self.jobs]
                multi_cmds = json.loads(self.commands)

                gmx = GMX(gmx_bin=Config.GMX_BIN)
                commands_list = gmx.generate_gpu_multidir_cmds(multi_dirs, multi_cmds,
                                                               n_parallel=Config.GMX_MULTI_NJOB,
                                                               n_gpu=jobmanager.ngpu,
                                                               n_omp=Config.GMX_MULTI_NOMP)
                for i, commands in enumerate(commands_list):
                    # instead of run directly, we add a record to pbs_job
                    sh = os.path.join(self.dir, '_job.run-%i.sh' % i)
                    pbs_name = '%s-run-%i' % (self.name, i)

                    if Config.GMX_MULTI_NOMP == None:
                        n_tasks = None
                    else:
                        n_tasks = Config.GMX_MULTI_NJOB * Config.GMX_MULTI_NOMP * \
                                  jobmanager.nprocs_request / jobmanager.nprocs

                    jobmanager.generate_sh(self.dir, commands, name=pbs_name, sh=sh, n_tasks=n_tasks)

                    pbs_job = PbsJob()
                    pbs_job.name = pbs_name
                    pbs_job.sh_file = sh
                    db.session.add(pbs_job)
                    db.session.flush()

                    # save pbs_job_id for jobs
                    for job in self.jobs[i * Config.GMX_MULTI_NJOB:(i + 1) * Config.GMX_MULTI_NJOB]:
                        job.pbs_job_id = pbs_job.id
                        job.status = Compute.Status.STARTED
                    db.session.commit()

                    # submit job, record if success or failed
                    pbs_job.submit()
                    time.sleep(0.2)

        except Exception as e:
            log.error('Run task failed %s %s' % (self, repr(e)))
            self.status = Compute.Status.FAILED
            db.session.commit()
            return False
        else:
            self.status = Compute.Status.STARTED
            db.session.commit()
            return True

    @property
    def ready_to_extend(self) -> bool:
        if not (self.stage == Compute.Stage.RUNNING and self.status == Compute.Status.STARTED):
            return False

        # do not extend the task if some job is running
        for job in self.jobs:
            if job.status == Compute.Status.STARTED:
                return False

        # do not extend the task if no job need extend
        for job in self.jobs:
            if job.need_extend:
                return True
        else:
            return False

    def extend(self, ignore_pbs_limit=False) -> bool:
        log.info('Extend %s' % self)

        if not self.ready_to_extend:
            log.warning('Not ready to extend %s' % self)
            return False

        jobs_extend = []
        for job in self.jobs:
            if job.need_extend:
                jobs_extend.append(job)
        if len(jobs_extend) == 0:
            log.warning('No job need to extend %s' % self)
            return False

        if self.cycle >= Config.EXTEND_CYCLE_LIMIT:
            log.warning('Will not extend %s EXTEND_CYCLE_LIMIT reached' % self)
            return False

        n_pbs_extend = len(jobs_extend)
        if Config.EXTEND_GMX_MULTI:
            n_pbs_extend = math.ceil(n_pbs_extend / Config.EXTEND_GMX_MULTI_NJOB)

        if not ignore_pbs_limit and jobmanager.n_running_jobs + n_pbs_extend > Config.PBS_NJOB_LIMIT:
            log.warning('PBS_NJOB_LIMIT reached')
            return False

        self.cycle += 1
        db.session.commit()

        try:
            if not Config.EXTEND_GMX_MULTI:
                for job in jobs_extend:
                    job.extend()
                    time.sleep(0.2)
            else:
                multi_dirs = []
                multi_cmds = []
                simulation = init_simulation(self.procedure, extend=True)
                for job in jobs_extend:
                    os.chdir(job.dir)
                    multi_dirs.append(job.dir)
                    multi_cmds = simulation.extend(jobname='%s-%i' % (job.name, job.cycle + 1),
                                                   sh='_job.extend-%i.sh' % (job.cycle + 1))

                commands_list = simulation.gmx.generate_gpu_multidir_cmds(multi_dirs, multi_cmds,
                                                                          n_parallel=Config.EXTEND_GMX_MULTI_NJOB,
                                                                          n_gpu=jobmanager.ngpu,
                                                                          n_procs=jobmanager.nprocs)

                os.chdir(self.dir)
                for i, commands in enumerate(commands_list):
                    sh = os.path.join(self.dir, '_job.extend-%i-%i.sh' % (self.cycle, i))
                    pbs_name = '%s-extend-%i-%i' % (self.name, self.cycle, i)

                    jm_extend.generate_sh(self.dir, commands, name=pbs_name, sh=sh)

                    # instead of run directly, we add a record to pbs_job
                    pbs_job = PbsJob()
                    pbs_job.name = pbs_name
                    pbs_job.sh_file = sh
                    db.session.add(pbs_job)
                    db.session.flush()

                    # save pbs_job_id for jobs
                    for job in jobs_extend[i * Config.EXTEND_GMX_MULTI_NJOB:(i + 1) * Config.EXTEND_GMX_MULTI_NJOB]:
                        job.pbs_job_id = pbs_job.id
                        job.cycle += 1
                        job.status = Compute.Status.STARTED
                    db.session.commit()

                    # submit job, record if success or failed
                    pbs_job.submit()
                    time.sleep(0.2)

        except Exception as e:
            log.error('Extend task failed %s %s' % (self, repr(e)))
            return False
        else:
            return True

    def check_finished(self):
        """
        check if all jobs in this tasks are finished
        if finished, analyze the job
        """
        log.info('Check status %s' % self)

        for job in self.jobs:
            try:
                job.check_finished()
            except Exception as e:
                log.error('Check job status failed %s %s' % (job, repr(e)))

            if job.status == Compute.Status.DONE:
                try:
                    log.info('Analyze %s' % job)
                    job.analyze()
                except Exception as e:
                    log.error('Analyze failed %s %s' % (job, repr(e)))

        # Set status as DONE only if all jobs are converged
        for job in self.jobs:
            if not job.converged:
                break
        else:
            self.status = Compute.Status.DONE
            db.session.commit()

    def check_finished_multiprocessing(self):
        """
        check if all jobs in this tasks are finished
        if finished, analyze the job
        """
        from multiprocessing import Pool

        log.info('Check status %s' % self)
        jobs_to_analyze = []
        for job in self.jobs:
            try:
                job.check_finished()
            except Exception as e:
                log.error('Check job status failed %s %s' % (job, repr(e)))

            if job.status == Compute.Status.DONE:
                log.info('Analyze %s' % job)
                jobs_to_analyze.append(job)

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
                if exception != None:
                    log.error('Analyze failed %s %s' % (job, exception))
                for k, v in return_dict.items():
                    setattr(job, k, v)

        db.session.commit()

        # Set status as DONE only if all jobs are converged
        for job in self.jobs:
            if not job.converged:
                break
        else:
            self.status = Compute.Status.DONE
            db.session.commit()

    def reset(self):
        for job in self.jobs:
            if not job.pbs_job_done:
                try:
                    jobmanager.kill_job(job.pbs_job.name)
                except Exception as e:
                    log.warning('Kill job %s Cannot kill PBS job: %s' % (job, repr(e)))
            db.session.delete(job)
        try:
            shutil.rmtree(self.dir)
        except:
            log.warning('Remove task %s Cannot remove folder: %s' % (self, self.dir))

        self.n_mol_list = None
        self.cycle = 0
        self.stage = Compute.Stage.SUBMITTED
        self.status = Compute.Status.DONE
        self.commands = None
        self.post_result = None
        db.session.commit()

    def remove(self):
        self.reset()
        db.session.delete(self)
        db.session.commit()

    def post_process(self, force=False):
        if not force:
            if not (self.stage == Compute.Stage.RUNNING and self.status == Compute.Status.DONE):
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

        if post_result != None:
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
    def pbs_job_done(self) -> bool:
        if self.pbs_job_id is None:
            return False
        elif not self.pbs_job.submitted:
            return False
        else:
            return not jobmanager.is_running(self.pbs_job.name)

    @property
    def need_extend(self) -> bool:
        return self.status == Compute.Status.ANALYZED and not self.converged

    def prepare(self):
        prior_job = self.prior_job
        if prior_job is None:
            prior_job_dir = None
        else:
            if not prior_job.converged:
                log.warning('Prepare job %s Prior job not converged' % repr(self))
                prior_job_dir = None
            else:
                prior_job_dir = prior_job.dir

        cd_or_create_and_cd(self.dir)
        simulation = init_simulation(self.task.procedure)
        commands = simulation.prepare(model_dir='../build', T=self.t, P=self.p, jobname=self.name,
                                      prior_job_dir=prior_job_dir, drde=True)  # Temperature dependent parameters

        # all jobs under one task share same commands, save this for GMX -multidir simulation
        if self.task.commands is None:
            self.task.commands = json.dumps(commands)
            db.session.commit()

    def run(self):
        try:
            os.chdir(self.dir)
        except:
            raise Exception('Should prepare job first')

        # instead of run directly, we add a record to pbs_job
        pbs_job = PbsJob()
        pbs_job.name = self.name
        pbs_job.sh_file = os.path.join(self.dir, jobmanager.sh)
        db.session.add(pbs_job)
        db.session.flush()

        self.pbs_job_id = pbs_job.id
        self.status = Compute.Status.STARTED
        db.session.commit()

        # submit job, record if success or failed
        pbs_job.submit()

    def extend(self):
        if not self.need_extend:
            log.warning('Will not extend %s Status: %s' % (self, Compute.Status.text[self.status]))
            return

        if self.cycle >= Config.EXTEND_CYCLE_LIMIT:
            log.warning('Will not extend %s EXTEND_CYCLE_LIMIT reached' % self)
            return

        os.chdir(self.dir)

        pbs_name = '%s-%i' % (self.name, self.cycle + 1)
        sh = os.path.join(self.dir, '_job.extend-%i.sh' % (self.cycle + 1))

        # instead of run directly, we add a record to pbs_job

        simulation = init_simulation(self.task.procedure, extend=True)
        simulation.extend(jobname=pbs_name, sh=sh)

        pbs_job = PbsJob()
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
        if not self.pbs_job_done:
            return False

        os.chdir(self.dir)

        simulation = init_simulation(self.task.procedure)
        if simulation.check_finished():
            self.status = Compute.Status.DONE
        else:
            self.status = Compute.Status.FAILED
            log.error('Job failed %s' % self)
        db.session.commit()
        return True

    def analyze(self, rerun=False, **kwargs):
        if self.status == Compute.Status.ANALYZED and not rerun:
            log.warning('Will not analyze %s Already analyzed' % self)
            return
        if self.status == Compute.Status.STARTED:
            log.warning('Will not analyze %s Job is still running' % self)
            return
        if self.status == Compute.Status.FAILED:
            log.warning('Will not analyze %s Job failed' % self)
            return

        os.chdir(self.dir)

        simulation = init_simulation(self.task.procedure)
        try:
            result = simulation.analyze(**kwargs)
        except:
            # One kine of fail -- simulation fail
            self.status = Compute.Status.FAILED
            db.session.commit()
            raise

        if result == None:
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

        # Clean intermediate files after analysis
        simulation.clean()
        db.session.commit()

    def analyze_multiprocessing(self, job_dir, job_procedure, **kwargs):
        _exception = None
        _status = self.status
        _converged = self.converged
        _result = self.result

        os.chdir(job_dir)
        simulation = init_simulation(job_procedure)
        try:
            result = simulation.analyze(**kwargs)
        except Exception as e:
            # One kine of fail -- simulation fail
            _status = Compute.Status.FAILED
            _exception = repr(e)
        else:
            if result == None:
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

        # Clean intermediate files after analysis
        simulation.clean()
        return {'exception': _exception,
                'status'   : _status,
                'converged': _converged,
                'result'   : _result,
                }

    def remove(self):
        if not self.pbs_job_done:
            try:
                jobmanager.kill_job(self.pbs_job.name)
            except Exception as e:
                log.warning('Kill job %s Cannot kill PBS job: %s' % (self, repr(e)))

        try:
            shutil.rmtree(self.dir)
        except:
            log.warning('Remove job %s Cannot remove folder: %s' % (self, self.dir))

        db.session.delete(self)
        db.session.commit()

    def remove_trj(self):
        os.chdir(self.dir)
        for f in os.listdir(os.getcwd()):
            if f.endswith('.trr') or f.endswith('.xtc') or f.startswith('.#'):
                try:
                    os.remove(f)
                except:
                    pass
