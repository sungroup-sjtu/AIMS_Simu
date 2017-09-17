import os
import socket
from collections import OrderedDict


class BaseConfig:
    CWD = os.path.dirname(os.path.abspath(__file__))
    DB_NAME = 'msdserver.sqlite'
    DB_PATH = os.path.join(CWD, DB_NAME)
    LOG = os.path.join(CWD, '_LOG_.txt')

    SQLALCHEMY_DATABASE_URI = 'sqlite:///%s' % DB_PATH
    SQLALCHEMY_COMMIT_ON_TEARDOWN = True
    SQLALCHEMY_TRACK_MODIFICATIONS = True

    MS_TOOLS_DIR = os.path.join(CWD, '../ms-tools')

    PBS_MANAGER = 'local'
    PBS_QUEUE_DICT = OrderedDict([(None, 2)])
    PBS_NJOB_LIMIT = 100
    PBS_ENV_CMD = ''

    GMX_MULTI = False  # do not perform gmx multidir simulation


class ClusterConfig(BaseConfig):
    WORK_DIR = '/share/workspace/msdserver/MSDServer/'
    DFF_ROOT = '/share/apps/dff/msdserver'
    PACKMOL_BIN = '/share/apps/tools/packmol'
    LMP_BIN = '/share/apps/lammps/lmp-stable'
    GMX_BIN = '/share/apps/gromacs/2016.3-msdserver/bin/gmx_gpu'

    PBS_MANAGER = 'torque'
    PBS_QUEUE_DICT = OrderedDict([('cpu', 8)])
    PBS_NJOB_LIMIT = 20


class TH2Config(BaseConfig):
    PBS_ENV_CMD = '''
source /BIGDATA/app/toolshs/cnmodule.sh
module purge
module load intel-compilers/mkl-15
module load gcc/5.3.0
'''

    WORK_DIR = '/BIGDATA/sjtu_hsun_1/MSDServer'
    DFF_ROOT = '/HOME/sjtu_hsun_1/apps/dff/7.3'
    PACKMOL_BIN = '/WORK/app/packmol/bin/packmol'
    GMX_BIN = '/HOME/sjtu_hsun_1/apps/gromacs/2016.3/bin/gmx_mpi'

    PBS_MANAGER = 'slurm'
    PBS_QUEUE_DICT = OrderedDict([('bigdata', 24)])
    PBS_NJOB_LIMIT = 64

    GMX_MULTI = True
    GMX_MULTI_NJOB = 8
    GMX_MULTI_NOMP = 6

    GMX_MULTI_EXTEND_NJOB = 8
    GMX_MULTI_EXTEND_SPECIAL_NJOB_NOMP = {1: 24, 2: 12, 3: 8, 5: 4, 6: 4}


class PIConfig(BaseConfig):
    PBS_ENV_CMD = '''
source /usr/share/Modules/init/bash
unset MODULEPATH
module use /lustre/usr/modulefiles/pi
module purge
module load icc/16.0 impi/5.1 mkl/11.3
export I_MPI_PMI_LIBRARY=/usr/lib64/libpmi.so
export I_MPI_FABRICS=shm:dapl
'''


class MacConfig(BaseConfig):
    WORK_DIR = '/tmp/MSDServer/'
    DFF_ROOT = '/Users/zheng/Projects/DFF/Developing'
    PACKMOL_BIN = '/Users/zheng/Projects/DFF/Developing/bin64m/Packmol/packmol.exe'
    LMP_BIN = '/Users/zheng/Projects/DFF/Developing/bin64m/Lammps/lammps'
    GMX_BIN = '/opt/gromacs/2016.3/bin/gmx'


Config = ClusterConfig
hostname = socket.gethostname()
if hostname == 'cluster.hpc.org':
    Config = ClusterConfig
elif hostname.startswith('ln') or hostname.startswith('cn'):
    Config = TH2Config
elif hostname == 'z-Mac.local':
    Config = MacConfig
else:
    raise Exception('msd-server will not work on this machine')
