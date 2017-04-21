import itertools
import json
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
    elif procedure == 'npt-cp':
        return simulationEngine.NptCp(**kwargs)
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
        smiles_list_components: [[str]] = detail['smiles_list']
        states = detail.get('states')

        for n, smiles_list in enumerate(itertools.product(*smiles_list_components)):
            if states == None:
                t_min = detail.get('t_min')
                t_max = detail.get('t_max')
                t_interval = detail.get('t_interval')
                t_number = detail.get('t_number')
                p_min = detail.get('p_min')
                p_max = detail.get('p_max')
            else:
                state = states[n]
                t_min = state.get('t_min')
                t_max = state.get('t_max')
                t_interval = state.get('t_interval')
                t_number = state.get('t_number')
                p_min = state.get('p_min')
                p_max = state.get('p_max')

            for procedure in detail['procedures']:
                task = Task()
                task.compute_id = self.id
                task.n_components = len(smiles_list_components)
                task.smiles_list = json.dumps(smiles_list)
                task.procedure = procedure
                if procedure in Procedure.T_RELEVANT:
                    task.t_min = t_min
                    task.t_max = t_max
                    task.t_interval = t_interval
                    task.t_number = t_number
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

        text = {
            STARTED: 'Started',
            DONE: 'Done',
            FAILED: 'Failed'
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
    t_number = Column(Integer, nullable=True)
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

        if self.t_min == None or self.t_max == None:
            T_list = [None]
        else:
            T_list = get_T_list_from_range(self.t_min, self.t_max, self.t_interval, self.t_number)
        if self.p_min == None or self.p_max == None:
            P_list = [None]
        else:
            P_list = get_P_list_from_range(self.p_min, self.p_max)

        for t in T_list:
            for p in P_list:
                if Job.query.filter(
                        and_(Job.task_id == self.id,
                             Job.t == t,
                             Job.p == p)
                ).first() != None:
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

    @property
    def dir(self) -> str:
        return os.path.join(Config.WORK_DIR, self.name)


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
    converged = Column(Boolean, nullable=True)
    result = Column(Text, nullable=True)
    next_cycle_started = NotNullColumn(Boolean, default=False)

    task = db.relationship(Task)

    def __repr__(self):
        return '<Job: %s %i-%i-%i>' % (self.name, self.t or 0, self.p or 0, self.cycle)

    def prepare(self):
        try:
            os.chdir(self.task.dir)
        except:
            raise Exception('Should build simulation box first')

        cd_or_create_and_cd(self.dir)
        simulation = init_simulation(self.task.procedure)
        simulation.prepare(model_dir='../build', T=self.t, P=self.p, nproc=Config.NPROC_PER_JOB, jobname=self.name)

    def run(self):
        try:
            os.chdir(self.dir)
        except:
            raise Exception('Should prepare job first')

        simulation = init_simulation(self.task.procedure)
        simulation.run()

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
        if self.status != Compute.Status.DONE:
            raise Exception('Job is still running')

        try:
            os.chdir(self.dir)
        except:
            raise

        simulation = init_simulation(self.task.procedure)
        try:
            converged, result = simulation.analyze()
        except Exception as e:
            print('Analyze failed: %s %s' % (self, str(e)))
        else:
            self.converged = converged
            if converged:
                self.result = json.dumps(result)
            db.session.commit()

    def start_next_cycle(self):
        if self.converged == None or self.converged:
            print('New cycle can only be started if this cycle does not converge')
            return

        if self.next_cycle_started:
            print('Next cycle already started')
            return

        self.next_cycle_started = True

        job = Job()
        job.task_id = self.task_id
        job.t = self.t
        job.p = self.p
        job.cycle = self.cycle + 1
        db.session.add(job)
        db.session.commit()

        job.prepare()
        job.run()

    @property
    def dir(self) -> str:
        dir_name = '%i-%i-%i' % (self.t or 0, self.p or 0, self.cycle)
        return os.path.join(self.task.dir, dir_name)

    @property
    def is_running(self) -> bool:
        return jobmanager.get_info_from_name(self.name)

    @property
    def next_cycle(self):
        return Job.query.filter(
            and_(
                Job.task_id == self.task_id,
                Job.t == self.t,
                Job.p == self.p,
                Job.cycle == self.cycle + 1
            )
        ).first()

    @property
    def previous_cycle(self):
        return Job.query.filter(
            and_(
                Job.task_id == self.task_id,
                Job.t == self.t,
                Job.p == self.p,
                Job.cycle == self.cycle - 1
            )
        ).first()


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
