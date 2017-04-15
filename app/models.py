import os
from datetime import datetime
from functools import partial
from threading import Thread

from sqlalchemy import Column, ForeignKey, Integer, Text, String, DateTime

from config import Config
from mstools.utils import random_string, cd_or_create_and_cd
from . import db

NotNullColumn = partial(Column, nullable=False)
from mstools.simulation import GmxSimulation, LammpsSimulation



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
    p_min = Column(Integer, nullable=True)
    p_max = Column(Integer, nullable=True)
    task_name = NotNullColumn(String(200), default=random_string)
    status = NotNullColumn(Integer, default=Compute.Status.SUBMITTED)
    remark = Column(Text, nullable=True)

    compute = db.relationship(Compute)
    jobs = db.relationship('Job', lazy='dynamic')

    def __repr__(self):
        return '<Task: %s %s>' % (self.smiles, self.procedure)

    def init_simulation(self):
        if Config.SIMULATION_ENGINE == 'lammps':
            self.simulation = LammpsSimulation(packmol_bin=Config.PACKMOL_BIN, dff_root=Config.DFF_ROOT,
                                               lmp_bin=Config.LMP_BIN, procedure=self.procedure)
        elif Config.SIMULATION_ENGINE == 'gmx':
            self.simulation = GmxSimulation(packmol_bin=Config.PACKMOL_BIN, dff_root=Config.DFF_ROOT,
                                            gmx_bin=Config.GMX_BIN, procedure=self.procedure)
        else:
            raise Exception('Simulation engine not supported')

    def build(self):
        cd_or_create_and_cd(self.base_dir)
        cd_or_create_and_cd('build')
        self.init_simulation()
        self.simulation.build(self.smiles, minimize=True)
        self.simulation.prepare()

    def run(self):
        try:
            os.chdir(self.base_dir)
        except:
            raise Exception('Should build simulation box first')
        else:
            self.init_simulation()
            self.simulation.run()

    def analyze(self):
        try:
            os.chdir(self.base_dir)
        except:
            raise Exception('Should build simulation box first')
        else:
            self.simulation.analyze()

    @property
    def base_dir(self) -> str:
        return os.path.join(Config.WORK_DIR, self.task_name)


class Job(db.Model):
    __tablename__ = 'job'
    id = NotNullColumn(Integer, primary_key=True)
    task_id = NotNullColumn(Integer, ForeignKey(Task.id))
    t = Column(Integer, nullable=True)
    p = Column(Integer, nullable=True)
    time = NotNullColumn(DateTime, default=datetime.now)
    status = NotNullColumn(Integer, default=Compute.Status.SUBMITTED)

    task = db.relationship(Task)

    def __repr__(self):
        return '<Job: %s %s>' % (self.smiles, self.procedure)

    @property
    def job_name(self):
        return self.task.task_name + '_' + self.t + '_' + self.p


class JobThread(Thread):
    def __init__(self, compute_id: int):
        super().__init__()
        self.compute_id = compute_id

    def run(self):
        compute = Compute.query.get(self.compute_id)
        Job = Task
        for job in Job.query.filter(Job.compute_id == self.compute_id):
            job.status = Compute.Status.PREPARING
            db.session.commit()
            job.build()
            job.status = Compute.Status.RUNNING
            db.session.commit()
            job.run_local()
            job.status = Compute.Status.DONE
            db.session.commit()
