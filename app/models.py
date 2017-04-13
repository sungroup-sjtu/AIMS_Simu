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


class ComputeProcedure:
    NPT = 'npt'
    NPT_SOLVATION = 'npt-solvation'
    NPT_HOV_IONIC_LIQUID = 'npt-hov-ionic-liquid'
    NVT = 'nvt'
    NVT_MSD = 'nvt-msd'
    NVT_VISCOSITY = 'nvt-viscosity'
    NVT_SLAB = 'nvt-slab'
    NVT_BINARY_SLAB = 'nvt-binary-slab'
    choices = [NPT, NPT_SOLVATION, NPT_HOV_IONIC_LIQUID, NVT, NVT_MSD, NVT_VISCOSITY, NVT_SLAB, NVT_BINARY_SLAB]
    T_RELEVANT = choices
    P_RELEVANT = [NPT, NPT_SOLVATION, NVT, NVT_MSD, NVT_VISCOSITY]


class Compute(db.Model):
    __tablename__ = 'compute'
    id = NotNullColumn(Integer, primary_key=True)
    web_id = NotNullColumn(Integer)
    web_user_id = NotNullColumn(Integer)
    web_ip = NotNullColumn(String(200))
    n_components = NotNullColumn(Integer)
    time = NotNullColumn(DateTime, default=datetime.now)
    json = NotNullColumn(Text)

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

    @property
    def jobs(self):
        if self.n_components == 1:
            Job = JobUnary
        elif self.n_components == 2:
            Job = JobBinary
        return Job.query.filter(Job.compute_id == self.id)


class JobUnary(db.Model):
    __tablename__ = 'job_unray'
    id = NotNullColumn(Integer, primary_key=True)
    compute_id = NotNullColumn(Integer, ForeignKey(Compute.id))
    smiles = NotNullColumn(Text)
    procedure = NotNullColumn(String(200))
    t = Column(Integer, nullable=True)
    p = Column(Integer, nullable=True)
    time = NotNullColumn(DateTime, default=datetime.now)
    job_name = NotNullColumn(String(200), default=random_string)
    status = NotNullColumn(Integer, default=Compute.Status.SUBMITTED)
    remark = Column(Text, nullable=True)
    n_mol = Column(Integer, nullable=True)

    compute = db.relationship(Compute)

    def __repr__(self):
        return '<Job: %s %s>' % (self.smiles, self.procedure)

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
        return os.path.join(Config.WORK_DIR, self.job_name)


class JobBinary(db.Model):
    __tablename__ = 'job_binary'
    id = NotNullColumn(Integer, primary_key=True)
    compute_id = NotNullColumn(Integer, ForeignKey(Compute.id))
    smiles = NotNullColumn(Text)
    smiles2 = NotNullColumn(Text)
    procedure = NotNullColumn(String(200))
    t = NotNullColumn(Integer)
    p = NotNullColumn(Integer)
    time = NotNullColumn(DateTime, default=datetime.now)
    job_name = Column(String(200), nullable=True)
    status = NotNullColumn(Integer, default=Compute.Status.SUBMITTED)
    remark = Column(Text, nullable=True)

    compute = db.relationship(Compute)

    def build(self):
        pass

    def prepare(self):
        pass

    def run(self):
        pass

    def analyze(self):
        pass


class JobThread(Thread):
    def __init__(self, compute_id: int):
        super().__init__()
        self.compute_id = compute_id

    def run(self):
        compute = Compute.query.get(self.compute_id)
        if compute.n_components == 1:
            Job = JobUnary
        elif compute.n_components == 2:
            Job = JobBinary
        for job in Job.query.filter(Job.compute_id == self.compute_id):
            job.status = Compute.Status.PREPARING
            db.session.commit()
            job.build()
            job.status = Compute.Status.RUNNING
            db.session.commit()
            job.run_local()
            job.status = Compute.Status.DONE
            db.session.commit()
