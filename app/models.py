import itertools
import json
import shutil
import time
from datetime import datetime
from functools import partial
from threading import Thread

from sqlalchemy import Column, ForeignKey, Integer, Text, String, Boolean, DateTime, and_

from config import Config
from mstools.simulation.procedure import Procedure
from mstools.utils import *
from . import db, jobmanager

NotNullColumn = partial(Column, nullable=False)


def init_simulation(procedure):
    if Config.SIMULATION_ENGINE == 'gmx':
        from mstools.simulation import gmx as simulationEngine
        kwargs = {'packmol_bin': Config.PACKMOL_BIN, 'dff_root': Config.DFF_ROOT,
                  'gmx_bin': Config.GMX_BIN, 'jobmanager': jobmanager}
    elif Config.SIMULATION_ENGINE == 'lammps':
        from mstools.simulation import lammps as simulationEngine
        kwargs = {'packmol_bin': Config.PACKMOL_BIN, 'dff_root': Config.DFF_ROOT,
                  'lmp_bin': Config.LMP_BIN, 'jobmanager': jobmanager}
    else:
        raise Exception('Simulation engine not supported')

    if procedure == 'npt':
        return simulationEngine.Npt(**kwargs)
    elif procedure == 'nvt':
        return simulationEngine.Nvt(**kwargs)
    elif procedure == 'nvt-slab':
        return simulationEngine.NvtSlab(**kwargs)
    elif procedure == 'nvt-vacuum':
        return simulationEngine.NvtVacuum(**kwargs)
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

    tasks = db.relationship('Task', lazy='dynamic')

    def create_tasks(self):
        detail = json.loads(self.json)
        procedures = detail['procedures']
        smiles_list_components: [[str]] = detail['smiles_list']
        states = detail.get('states')

        # check for prior
        for p in detail['procedures']:
            prior = Procedure.prior.get(p)
            if prior != None and prior not in procedures:
                procedures = [prior] + procedures  # Put prerequisite at first

        for n, smiles_list in enumerate(itertools.product(*smiles_list_components)):
            if states == None:
                state = detail.get('state')
            else:
                state = states[n]

            t_min = state.get('t_min')
            t_max = state.get('t_max')
            t_interval = state.get('t_interval')
            p_min = state.get('p_min')
            p_max = state.get('p_max')

            for procedure in procedures:
                task = Task()
                task.compute_id = self.id
                task.n_components = len(smiles_list_components)
                task.smiles_list = json.dumps(smiles_list)
                task.procedure = procedure
                if procedure in Procedure.T_RELEVANT:
                    task.t_min = t_min
                    task.t_max = t_max
                    task.t_interval = t_interval
                    if t_min == None or t_max == None:
                        raise Exception('Invalid temperature')
                else:
                    task.t_min = None
                    task.t_max = None
                if procedure in Procedure.P_RELEVANT:
                    task.p_min = p_min
                    task.p_max = p_max
                    if p_min == None or p_max == None:
                        raise Exception('Invalid pressure')
                else:
                    task.p_min = None
                    task.p_max = None
                db.session.add(task)

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
        return '<Task: %s %s>' % (self.smiles_list, self.procedure)

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
            simulation.build(minimize=True)
            self.n_mol_list = json.dumps(simulation.n_mol_list)
        except:
            self.status = Compute.Status.FAILED
            db.session.commit()
            raise
        else:
            try:
                self.prepare_jobs()
            except:
                self.status = Compute.Status.FAILED
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
            T_list = get_T_list_from_range(self.t_min, self.t_max, self.t_interval)
        if self.p_min is None or self.p_max is None:
            P_list = [None]
        else:
            P_list = get_P_list_from_range(self.p_min, self.p_max)

        for t in T_list:
            for p in P_list:
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
                db.session.add(job)

        try:
            db.session.commit()
        except:
            db.session.rollback()
            raise Exception('Prepare jobs failed')
        else:
            for job in self.jobs:
                job.prepare()

    def run(self, sleep=0.1):
        if self.stage == Compute.Stage.BUILDING and self.status == Compute.Status.DONE:
            self.stage = Compute.Stage.RUNNING
            self.status = Compute.Status.STARTED
            db.session.commit()
            for job in self.jobs:
                job.run()
                time.sleep(sleep)
        else:
            raise Exception('Should build first')

    def check_finished(self):
        finished = True
        failed = False
        for job in self.jobs:
            if not job.check_finished():
                finished = False
            elif job.status == Compute.Status.FAILED:
                failed = True
            elif not job.converged:
                finished = False

        if finished:
            if failed:
                self.status = Compute.Status.FAILED
            else:
                self.status = Compute.Status.DONE
        db.session.commit()

    def remove(self):
        for job in self.jobs:
            job.remove()
        try:
            shutil.rmtree(self.dir)
        except:
            print('Cannot remove folder: %s' % self.dir)
        else:
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

    task = db.relationship(Task)

    def __repr__(self):
        return '<Job: %s %i-%i-%i>' % (self.name, self.t or 0, self.p or 0, self.cycle)

    @property
    def dir(self) -> str:
        dir_name = '%i-%i-%i' % (self.t or 0, self.p or 0, self.cycle)
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
        return jobmanager.get_info(self.name)

    def prepare(self):
        prior_job: Job = self.prior_job
        if prior_job != None:
            if not prior_job.converged:
                print('Prior job has not finished, this job will not be prepared now')
                return
            else:
                prior_job_dir = prior_job.dir
                prior_job_result = prior_job.result
        else:
            prior_job_dir = None
            prior_job_result = None

        try:
            os.chdir(self.task.dir)
        except:
            raise Exception('Should build simulation box first')

        cd_or_create_and_cd(self.dir)
        simulation = init_simulation(self.task.procedure)
        simulation.prepare(model_dir='../build', T=self.t, P=self.p, jobname=self.name,
                           prior_job_dir=prior_job_dir, prior_job_result=prior_job_result)

    def run(self):
        try:
            os.chdir(self.dir)
        except:
            raise Exception('Should prepare job first')

        simulation = init_simulation(self.task.procedure)
        simulation.run()

    def extend(self):
        if self.status == Compute.Status.STARTED:
            print('Job is still running')
            return
        if self.status == Compute.Status.FAILED:
            print('Job failed, will not extend')
            return
        if self.converged:
            print('Already converged, no need to extend')
            return

        try:
            os.chdir(self.dir)
        except:
            raise

        simulation = init_simulation(self.task.procedure)
        simulation.extend(jobname=self.name)
        simulation.run()

        self.cycle += 1
        self.status = Compute.Status.STARTED
        db.session.commit()

    def check_finished(self) -> bool:
        if self.status in (Compute.Status.DONE, Compute.Status.FAILED):
            return True
        if self.is_running:
            return False

        try:
            os.chdir(self.dir)
        except:
            raise Exception('Should prepare job first')

        simulation = init_simulation(self.task.procedure)
        if simulation.check_finished():
            self.status = Compute.Status.DONE
        else:
            self.status = Compute.Status.FAILED
        db.session.commit()
        return True

    def analyze(self):
        if self.converged:
            print('Already converged')
            return
        if self.status == Compute.Status.STARTED:
            print('Job is still running')
            return
        if self.status == Compute.Status.FAILED:
            print('Job failed, will not perform analyze')
            return

        try:
            os.chdir(self.dir)
        except:
            raise

        simulation = init_simulation(self.task.procedure)
        dirs = [self.dir]
        try:
            result = simulation.analyze(dirs)
        except Exception as e:
            print('Analyze failed: %s %s' % (self, str(e)))
            self.status = Compute.Status.FAILED
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
            jobmanager.kill_job(self.name)
        except Exception as e:
            print(str(e))

        try:
            shutil.rmtree(self.dir)
        except:
            print('Cannot remove folder: %s' % self.dir)
        else:
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
