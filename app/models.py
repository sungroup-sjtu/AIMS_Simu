import itertools
import json
from datetime import datetime
from functools import partial
from threading import Thread

from sqlalchemy import Column, ForeignKey, Integer, Text, String, DateTime, and_

from config import Config
from mstools.simulation.procedure import Procedure
from mstools.utils import *
from . import db, simulation

NotNullColumn = partial(Column, nullable=False)


class Compute(db.Model):
    __tablename__ = 'compute'
    id = NotNullColumn(Integer, primary_key=True)
    web_id = NotNullColumn(Integer)
    web_user_id = NotNullColumn(Integer)
    web_ip = NotNullColumn(String(200))
    time = NotNullColumn(DateTime, default=datetime.now)
    json = NotNullColumn(Text)

    tasks = db.relationship('Task', lazy='dynamic')

    class Status:
        SUBMITTED = 0
        PREPARING = 1
        RUNNING = 2
        DONE = 9
        FAILED = -1

        text = {
            SUBMITTED: 'Submitted',
            PREPARING: 'Building...',
            RUNNING: 'Running...',
            DONE: 'Done',
            FAILED: 'Failed'
        }

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
    status = NotNullColumn(Integer, default=Compute.Status.SUBMITTED)
    remark = Column(Text, nullable=True)

    compute = db.relationship(Compute)
    jobs = db.relationship('Job', lazy='dynamic')

    def __repr__(self):
        return '<Task: %s %s>' % (self.smiles_list, self.procedure)

    def build(self):
        cd_or_create_and_cd(self.dir)
        cd_or_create_and_cd('build')
        simulation.set_procedure(self.procedure)
        simulation.build(json.loads(self.smiles_list), minimize=True)
        self.prepare_jobs()

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

    def run(self):
        for job in self.jobs:
            job.run()

    def analyze(self):
        for job in self.jobs:
            job.analyze()

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
    status = NotNullColumn(Integer, default=Compute.Status.SUBMITTED)
    name = NotNullColumn(String(200), default=random_string)

    task = db.relationship(Task)

    def __repr__(self):
        return '<Job: %s %i-%i>' % (self.name, self.t or 0, self.p or 0)

    def prepare(self):
        try:
            os.chdir(self.task.dir)
        except:
            raise Exception('Should build simulation box first')

        cd_or_create_and_cd(self.dir)
        simulation.set_procedure(self.task.procedure)
        simulation.prepare(model_dir='../build', T=self.t, P=self.p, nproc=Config.NPROC_PER_JOB, jobname=self.name)

    def run(self):
        try:
            os.chdir(self.dir)
        except:
            raise Exception('Should prepare job first')

        simulation.run()

    def analyze(self):
        pass

    @property
    def dir(self) -> str:
        dir_name = '%i-%i' % (self.t or 0, self.p or 0)
        return os.path.join(self.task.dir, dir_name)


class TaskThread(Thread):
    def __init__(self, app, compute_id: int):
        super().__init__()
        self.app = app
        self.compute_id = compute_id

    def run(self):
        with self.app.app_context():
            for task in Task.query.filter(Task.compute_id == self.compute_id):
                task.status = Compute.Status.PREPARING
                db.session.commit()
                task.build()
                task.status = Compute.Status.RUNNING
                db.session.commit()
                task.run()
