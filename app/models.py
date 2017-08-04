import sys
import re
import json
import shutil
import time
import warnings
from datetime import datetime
from functools import partial
from threading import Thread

from sqlalchemy import Column, ForeignKey, Integer, Text, String, Boolean, DateTime, and_

from . import db

from config import Config

sys.path.append(Config.MS_TOOLS_DIR)
from mstools.simulation.procedure import Procedure
from mstools.utils import *

NotNullColumn = partial(Column, nullable=False)

from mstools.jobmanager import *
from mstools.wrapper import GMX

if Config.PBS_MANAGER == 'local':
    jobmanager = Local(queue_dict=Config.PBS_QUEUE_DICT, env_cmd=Config.PBS_ENV_CMD)
elif Config.PBS_MANAGER == 'torque':
    jobmanager = Torque(queue_dict=Config.PBS_QUEUE_DICT, env_cmd=Config.PBS_ENV_CMD)
elif Config.PBS_MANAGER == 'slurm':
    jobmanager = Slurm(queue_dict=Config.PBS_QUEUE_DICT, env_cmd=Config.PBS_ENV_CMD)
else:
    raise Exception('Job manager not supported')


def init_simulation(procedure):
    from mstools.simulation import gmx as simulationEngine
    kwargs = {'packmol_bin': Config.PACKMOL_BIN, 'dff_root': Config.DFF_ROOT,
              'gmx_bin': Config.GMX_BIN, 'jobmanager': jobmanager}

    if procedure == 'npt':
        return simulationEngine.Npt(**kwargs)
    elif procedure == 'nvt':
        return simulationEngine.Nvt(**kwargs)
    elif procedure == 'nvt-slab':
        return simulationEngine.NvtSlab(**kwargs)
    else:
        raise Exception('Unknown simulation procedure')


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

    def create_tasks(self):
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

                prior = Procedure.prior.get(p)
                if prior is not None and prior not in procedures:
                    procedures = [prior] + procedures  # Put prerequisite at first

            all_smiles_list = []
            for combination in combinations:
                # check for smiles and name
                smiles_list = combination['smiles']
                name_list = combination.get('names')

                if smiles_list in all_smiles_list:
                    raise Exception('Duplicated smiles combinations')
                all_smiles_list.append(smiles_list)

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
                    if Task.query.filter(
                            and_(Task.compute_id == self.id,
                                 Task.smiles_list == json.dumps(smiles_list),
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
            raise Exception('Cannot create tasks: ' + repr(e))

    class Stage:
        SUBMITTED = 0
        BUILDING = 1
        RUNNING = 2

        text = {
            SUBMITTED: 'Submitted',
            BUILDING: 'Building...',
            RUNNING: 'Running...',
        }

    class Status:
        STARTED = 1
        DONE = 9
        FAILED = -1
        ANALYZED = 10

        text = {
            STARTED: 'Started',
            DONE: 'Done',
            FAILED: 'Failed',
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
    stage = NotNullColumn(Integer, default=Compute.Stage.SUBMITTED)
    status = NotNullColumn(Integer, default=Compute.Status.DONE)
    remark = Column(Text, nullable=True)

    compute = db.relationship(Compute)
    jobs = db.relationship('Job', lazy='dynamic')

    def __repr__(self):
        return '<Task: %s %s %s>' % (self.smiles_list, self.procedure, self.name)

    @property
    def dir(self) -> str:
        return os.path.join(Config.WORK_DIR, self.name)

    @property
    def prior_task(self):
        prior_procedure = Procedure.prior.get(self.procedure)
        if prior_procedure == None:
            return None

        return Task.query.filter(
            and_(Task.smiles_list == self.smiles_list,
                 Task.procedure == prior_procedure,
                 Task.t_min == self.t_min,
                 Task.t_max == self.t_max,
                 Task.t_interval == self.t_interval,
                 Task.p_min == self.p_max,
                 Task.p_max == self.p_max
                 )
        ).first()

    @property
    def n_pbs_jobs(self) -> int:
        n = len(self.jobs.all())
        if Config.GMX_MULTIDIR:
            n = math.ceil(n / Config.GMX_MULTIDIR_NJOB)
        return n

    def build(self):
        # TODO optimize logic
        cd_or_create_and_cd(self.dir)
        cd_or_create_and_cd('build')

        self.stage = Compute.Stage.BUILDING
        self.status = Compute.Status.STARTED
        db.session.commit()

        simulation = init_simulation(self.procedure)
        try:
            simulation.set_system(json.loads(self.smiles_list), n_atoms=3000)
            simulation.build()
            self.n_mol_list = json.dumps(simulation.n_mol_list)
        except Exception as e:
            self.status = Compute.Status.FAILED
            self.remark = repr(e)
            db.session.commit()
            raise
        else:
            try:
                self.prepare_jobs()
            except Exception as e:
                self.status = Compute.Status.FAILED
                self.remark = repr(e)
                db.session.commit()
                raise
            else:
                self.status = Compute.Status.DONE
                db.session.commit()

    def prepare_jobs(self):
        if not os.path.exists(os.path.join(self.dir, 'build')):
            raise Exception('Should build simulation box first')

        if self.t_min is None or self.t_max is None:
            T_list = [None]
        else:
            T_list = get_T_list_from_range(self.t_min, self.t_max, interval=self.t_interval)
        if self.p_min is None or self.p_max is None:
            P_list = [None]
        else:
            P_list = get_P_list_from_range(self.p_min, self.p_max, multiple=(2, 5))
            # remove P in the range of (1,10] bar
            P_list = list(filter(lambda x: x == int(1E5) or x > int(1E6), P_list))

        for p in P_list:
            for t in T_list:
                if Job.query.filter(
                        and_(Job.task_id == self.id,
                             Job.t == t,
                             Job.p == p)
                ).first() is not None:
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

        try:
            for job in self.jobs:
                job.prepare()
        except:
            with warnings.catch_warnings(record=True):  # do not show warnings
                self.remove()
            raise

    def run(self, ignore_pbs_limit=False, sleep=0.2):
        if self.stage != Compute.Stage.BUILDING:
            raise Exception('Incorrect stage: %s' % Compute.Stage.text[self.stage])
        elif self.status != Compute.Status.DONE:
            raise Exception('Incorrect status: %s' % Compute.Status.text[self.status])

        if not ignore_pbs_limit and jobmanager.n_running_jobs + self.n_pbs_jobs >= Config.PBS_NJOB_LIMIT:
            raise Exception('PBS_NJOB_LIMIT reached, will not run job now')

        os.chdir(self.dir)

        self.stage = Compute.Stage.RUNNING
        self.status = Compute.Status.STARTED
        db.session.commit()

        if not Config.GMX_MULTIDIR:
            for job in self.jobs:
                job.run()
                time.sleep(sleep)
        else:
            multi_dirs = []
            multi_cmds = []
            for job in self.jobs:
                multi_dirs.append(job.dir)
                multi_cmds = json.loads(job.commands)
            gmx = GMX(gmx_bin=Config.GMX_BIN)
            commands_list = gmx.generate_gpu_multidir_cmds(multi_dirs, multi_cmds,
                                                           n_parallel=Config.GMX_MULTIDIR_NJOB,
                                                           n_gpu=0,
                                                           n_thread=Config.GMX_MULTIDIR_NTHREAD)
            for i, commands in enumerate(commands_list):
                sh = os.path.join(self.dir, '_job.multi-%i.sh' % i)
                pbs_name = '%s-%i' % (self.name, i)
                jobmanager.generate_sh(self.dir, commands, name=pbs_name, sh=sh,
                                       n_thread=Config.GMX_MULTIDIR_NTHREAD, exclusive=True)
                jobmanager.submit(sh)

                # save pbs_name for jobs
                for job in self.jobs[i * Config.GMX_MULTIDIR_NJOB:(i + 1) * Config.GMX_MULTIDIR_NJOB]:
                    job.pbs_name = pbs_name
                db.session.commit()
                time.sleep(sleep)

    def check_finished(self):
        '''
        set status as FAILED if one job failed
        set status as DONE only if all jobs are converged
        '''
        for job in self.jobs:
            try:
                job.check_finished()
            except Exception as e:
                warnings.warn('Check job status failed %s %s' % (job, repr(e)))

        for job in self.jobs:
            if job.status == Compute.Status.FAILED:
                self.status = Compute.Status.FAILED
                db.session.commit()
                warnings.warn('Task %s has failed job %s' % (self, job))
                return

        for job in self.jobs:
            if not job.converged:
                self.status = Compute.Status.STARTED
                db.session.commit()
                return

        self.status = Compute.Status.DONE
        db.session.commit()

    def remove(self):
        for job in self.jobs:
            job.remove()
        try:
            shutil.rmtree(self.dir)
        except:
            warnings.warn('Remove task %s Cannot remove folder: %s' % (self, self.dir))

        db.session.delete(self)
        db.session.commit()


class Job(db.Model):
    __tablename__ = 'job'
    id = NotNullColumn(Integer, primary_key=True)
    task_id = NotNullColumn(Integer, ForeignKey(Task.id))
    t = Column(Integer, nullable=True)
    p = Column(Integer, nullable=True)
    time = NotNullColumn(DateTime, default=datetime.now)
    name = NotNullColumn(String(200), default=random_string)
    cycle = NotNullColumn(Integer, default=1)
    status = NotNullColumn(Integer, default=Compute.Status.STARTED)
    converged = NotNullColumn(Boolean, default=False)
    result = Column(Text, nullable=True)
    commands = Column(Text, nullable=True)
    pbs_name = Column(String(200), nullable=True)

    task = db.relationship(Task)

    def __repr__(self):
        return '<Job: %s %i>' % (self.name, self.cycle)

    @property
    def dir(self) -> str:
        dir_name = '%i-%i' % (self.t or 0, self.p or 0)
        return os.path.join(self.task.dir, dir_name)

    @property
    def prior_job(self):
        prior_task = self.task.prior_task
        if prior_task == None:
            return None

        for job in prior_task.jobs:
            if job.t == self.t and job.p == self.p:
                return job

    @property
    def is_running(self) -> bool:
        return jobmanager.is_running(self.pbs_name)

    def prepare(self):
        prior_job = self.prior_job
        if prior_job is not None:
            if not prior_job.converged:
                warnings.warn('Prepare job %s Prior job not converged' % repr(self))
                return
            else:
                prior_job_dir = prior_job.dir
        else:
            prior_job_dir = None

        cd_or_create_and_cd(self.dir)
        simulation = init_simulation(self.task.procedure)
        commands = simulation.prepare(model_dir='../build', T=self.t, P=self.p, jobname=self.name,
                                      prior_job_dir=prior_job_dir, drde=True)  # Temperature dependent parameters
        self.commands = json.dumps(commands)
        db.session.commit()

    def run(self):
        try:
            os.chdir(self.dir)
        except:
            raise Exception('Should prepare job first')

        simulation = init_simulation(self.task.procedure)
        simulation.run()

        self.pbs_name = self.name
        db.session.commit()

    def extend(self, ignore_pbs_limit=False):
        if self.status == Compute.Status.STARTED:
            warnings.warn('Will not extend %s Job still running' % self)
            return
        if self.status == Compute.Status.FAILED:
            warnings.warn('Will not extend %s Job failed' % self)
            return
        if self.converged:
            warnings.warn('Will not extend %s Job converged' % self)
            return

        if not ignore_pbs_limit and jobmanager.n_running_jobs + 1 >= Config.PBS_NJOB_LIMIT:
            raise Exception('PBS_NJOB_LIMIT reached, will not extend job now')

        try:
            os.chdir(self.dir)
        except:
            raise

        pbs_name = '%s-%i' % (self.name, self.cycle)
        simulation = init_simulation(self.task.procedure)
        simulation.extend(jobname=pbs_name)
        simulation.run()

        self.cycle += 1
        self.status = Compute.Status.STARTED
        self.pbs_name = pbs_name
        db.session.commit()

    def check_finished(self) -> bool:
        if self.status in (Compute.Status.DONE, Compute.Status.FAILED, Compute.Status.ANALYZED):
            return True
        if self.is_running:
            return False

        simulation = init_simulation(self.task.procedure)
        if simulation.check_finished():
            self.status = Compute.Status.DONE
        else:
            self.status = Compute.Status.FAILED
        db.session.commit()
        return True

    def analyze(self, rerun=False):
        if self.status == Compute.Status.ANALYZED and not rerun:
            warnings.warn('Will not analyze %s Already analyzed' % self)
            return
        if self.status == Compute.Status.STARTED:
            warnings.warn('Will not analyze %s Job is still running' % self)
            return
        if self.status == Compute.Status.FAILED:
            warnings.warn('Will not analyze %s Job failed' % self)
            return

        try:
            os.chdir(self.dir)
        except:
            raise

        simulation = init_simulation(self.task.procedure)
        dirs = [self.dir]
        try:
            result = simulation.analyze(dirs)
        except:
            self.status = Compute.Status.FAILED
            db.session.commit()
            raise
        else:
            self.status = Compute.Status.ANALYZED
            if result is None:
                self.converged = False
                self.result = None
            else:
                self.converged = True
                self.result = json.dumps(result)
            db.session.commit()

    def remove(self):
        try:
            jobmanager.kill_job(self.pbs_name)
        except Exception as e:
            warnings.warn('Remove job %s Cannot kill PBS job: %s' % (self, repr(e)))

        try:
            shutil.rmtree(self.dir)
        except:
            warnings.warn('Remove job %s Cannot remove folder: %s' % (self, self.dir))

        db.session.delete(self)
        db.session.commit()


class TaskThread(Thread):
    def __init__(self, app, compute_id: int):
        super().__init__()
        self.app = app
        self.compute_id = compute_id

    def run(self):
        with self.app.app_context():
            for task in Task.query.filter(Task.compute_id == self.compute_id):
                task.build()
                task.run()
