import os
import socket


class BaseConfig:
    DB_NAME = 'msdserver.sqlite'
    DB_PATH = os.path.dirname(os.path.abspath(__file__)) + os.sep + DB_NAME

    SQLALCHEMY_DATABASE_URI = 'sqlite:///%s' % DB_PATH
    SQLALCHEMY_COMMIT_ON_TEARDOWN = True
    SQLALCHEMY_TRACK_MODIFICATIONS = True

    SIMULATION_ENGINE = 'gmx'
    # SIMULATION_ENGINE = 'lammps'

    JOB_MANAGER = 'local'
    NPROC_PER_JOB = 1

    # JOB_MANAGER = 'torque'
    # JOB_QUEUE = 'cpu'
    # NPROC_PER_JOB = 8


class ClusterConfig(BaseConfig):
    USER = 'msdserver'
    WORK_DIR = '/share/workspace/msdserver/MSDataServer/'
    DFF_ROOT = '/share/apps/dff/msdserver'
    PACKMOL_BIN = '/share/apps/tools/packmol'
    LMP_BIN = '/share/apps/lammps/lmp-stable'
    GMX_BIN = '/share/apps/gromacs/2016.2/bin/gmx'


class MacConfig(BaseConfig):
    USER = 'zheng'
    WORK_DIR = '/tmp/MSDataServer/'
    DFF_ROOT = '/Users/zheng/Projects/DFF7.2'
    PACKMOL_BIN = '/Users/zheng/Projects/DFF7.2/bin32m/Packmol/packmol.exe'
    LMP_BIN = '/Users/zheng/Projects/DFF7.2/bin32m/Lammps/lammps'
    GMX_BIN = '/usr/local/bin/gmx'


class XPSConfig(BaseConfig):
    USER = 'zgong'
    WORK_DIR = 'D:/Download/Temp'
    DFF_ROOT = 'D:/Projects/DFF7.2'
    PACKMOL_BIN = 'D:/Projects/DFF7.2/bin32w/Packmol/packmol.exe'
    LMP_BIN = 'D:/Projects/DFF7.2/bin32w/Lammps/lammps.exe'
    GMX_BIN = None


Config = ClusterConfig
if socket.gethostname() == 'z-Mac.local':
    Config = MacConfig
elif socket.gethostname() == 'z-XPS':
    Config = XPSConfig
else:
    raise Exception('MSData Server will not work on this machine')

