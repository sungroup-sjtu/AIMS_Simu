import os
import shutil
from datetime import datetime
from functools import partial

import math
from sqlalchemy import Column, ForeignKey, Integer, Text, String, DateTime

from config import Config
from mstools.tools import random_string, cd_or_create_and_cd
from . import db

NotNullColumn = partial(Column, nullable=False)
from mstools.wrapper import DFF, Packmol, Lammps
from mstools.simulation import Simulation


class Unit:
    K = 1
    Pa = 1
    kPa = 1000
    MPa = int(1E6)
    GPa = int(1E9)
    bar = int(1E5)
    atm = int(1.013E5)
    J = 1
    kJ = 1000
    kcal = 4184
    J_per_mol = 1
    kJ_per_mol = 1000
    kcal_per_mol = 4184

    text = {
        K: 'K',
        Pa: 'Pa',
        kPa: 'kPa',
        MPa: 'MPa',
        GPa: 'GPa',
        bar: 'bar',
        atm: 'atm',
        J: 'J',
        kJ: 'kJ',
        kcal: 'kcal',
        J_per_mol: 'J/mol',
        kJ_per_mol: 'kJ/mol',
        kcal_per_mol: 'kcal/mol'
    }


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

    def build(self):
        cd_or_create_and_cd(self.base_dir)
        cd_or_create_and_cd('build')

        simulation = Simulation(packmol_bin=Config.PACKMOL_BIN, dff_root=Config.DFF_ROOT, lmp_bin=Config.LAMMPS_BIN)
        simulation.build_lammps_box_from_smiles(self.smiles, 3000, 'init.data', 'em.lmp', minimize=True)

        print('Preparing simulation files...')
        os.chdir('..')
        shutil.copy('build/em.data', 'em.data')
        special, bond, angle, dihedral, improper = Lammps.get_intra_style_from_lmp('build/em.lmp')
        Lammps.prepare_lmp_from_template('t_npt.lmp', 'in.lmp', 'em.data', self.t, self.p / Unit.bar, int(1E3),
                                         self.n_mol, special, bond, angle, dihedral, improper)

    def run_local(self):
        try:
            os.chdir(self.base_dir)
        except:
            raise Exception('Should build simulation box first')

        if not (os.path.exists('in.lmp') and os.path.exists('em.data')):
            raise Exception('Should prepare simulation first')

        lammps = Lammps(Config.LAMMPS_BIN)
        print('Running NPT simulation...')
        lammps.run('in.lmp')

    def analyze(self):
        pass

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
