import os
import socket


class BaseConfig:
    CWD = os.path.dirname(os.path.abspath(__file__))
    DB = 'database/msdserver.sqlite'
    LOG = os.path.join(CWD, '_LOG_.txt')

    SQLALCHEMY_DATABASE_URI = 'sqlite:///%s?check_same_thread=False' % os.path.join(CWD, DB)
    SQLALCHEMY_BINDS = {
        'cv'  : 'sqlite:///%s?check_same_thread=False' % os.path.join(CWD, 'database/cv.sqlite'),
        'yaws': 'sqlite:///%s?check_same_thread=False' % os.path.join(CWD, 'database/yaws.sqlite'),
        'nist': 'sqlite:///%s?check_same_thread=False' % os.path.join(CWD, 'database/nist.sqlite'),
    }
    SQLALCHEMY_COMMIT_ON_TEARDOWN = True
    SQLALCHEMY_TRACK_MODIFICATIONS = True

    MS_TOOLS_DIR = os.path.join(CWD, '../ms-tools')

    PBS_MANAGER = 'local'
    PBS_QUEUE_LIST = [(2, 0)]
    PBS_NJOB_LIMIT = 10
    PBS_ENV_CMD = ''

    GMX_MDRUN = None
    EXTEND_GMX_MDRUN = None
    GMX_MULTI = False  # do not perform gmx multidir simulation

    EXTEND_CYCLE_LIMIT = 8

    DFF_TABLE = 'MGI'


class ClusterConfig(BaseConfig):
    PBS_ENV_CMD = '''
module purge
module load gcc openmpi gromacs/2016.5
'''

    WORK_DIR = '/share/workspace/gongzheng/_MSDServer/'
    DFF_ROOT = '/home/gongzheng/apps/DFF/Developing'
    PACKMOL_BIN = '/share/apps/tools/packmol'

    PBS_MANAGER = 'slurm'
    PBS_NJOB_LIMIT = 48
    PBS_QUEUE_LIST = [('gtx', 32, 2, 16)]  # partition, cpu(hyperthreading), gpu, cpu_request
    PBS_TIME_LIMIT = 1.0  # hour

    GMX_BIN = 'gmx_gpu'
    GMX_MULTI = True
    GMX_MULTI_NJOB = 8
    GMX_MULTI_NOMP = None  # Use only one node. Automatically determine the best number of threads.

    ### Extend
    EXTEND_PBS_QUEUE_LIST = [('fast', 24, 0, 12)]
    EXTEND_GMX_BIN = 'gmx_fast'
    EXTEND_GMX_MULTI = False  # Do not run -multidir simulation for Extend. So each job have same length
    EXTEND_GMX_MULTI_NJOB = 2  # Not used


class PIConfig(BaseConfig):
    PBS_ENV_CMD = '''
source /usr/share/Modules/init/bash
module purge
unset MODULEPATH
module use /lustre/usr/modulefiles/pi
module load icc/16.0 mkl/11.3 impi/5.1 tbb/4.3
export I_MPI_PMI_LIBRARY=/usr/lib64/libpmi.so
export I_MPI_FABRICS=shm:dapl
'''

    _base_dir = '/lustre/home/acct-nishsun/nishsun-1/'
    WORK_DIR = _base_dir + 'workspace/_MSDServer/'
    DFF_ROOT = _base_dir + 'apps/DFF/Developing'
    PACKMOL_BIN = _base_dir + 'apps/tools/packmol'

    PBS_MANAGER = 'slurm'
    PBS_SUBMIT_CMD = 'sbatch --reservation=cpu_nishsun'
    PBS_NJOB_LIMIT = 64
    PBS_QUEUE_LIST = [('cpu', 16, 0, 16)]  # partition, cpu(hyperthreading), gpu, cpu_request

    GMX_BIN = 'gmx_serial'
    GMX_MDRUN = 'gmx mdrun'
    GMX_MULTI = True
    GMX_MULTI_NJOB = 4
    GMX_MULTI_NOMP = 4

    ### Extend
    EXTEND_PBS_QUEUE_LIST = [('cpu', 16, 0, 16)]
    EXTEND_GMX_BIN = 'gmx_serial'
    EXTEND_GMX_MDRUN = 'gmx mdrun'
    EXTEND_GMX_MULTI = False  # Do not run -multidir simulation for Extend. So each job have same length
    EXTEND_GMX_MULTI_NJOB = 2  # Not used


class TH2Config(BaseConfig):
    PBS_ENV_CMD = '''
source /BIGDATA/app/toolshs/cnmodule.sh
module purge
module load intel-compilers/mkl-15
module load gcc/5.3.0
'''

    WORK_DIR = '/BIGDATA/sjtu_hsun_1/_MSDServer'
    DFF_ROOT = '/HOME/sjtu_hsun_1/apps/DFF/Developing'
    PACKMOL_BIN = '/WORK/app/packmol/bin/packmol'
    GMX_BIN = '/HOME/sjtu_hsun_1/apps/gromacs/2016.3/bin/gmx_mpi'

    PBS_MANAGER = 'slurm'
    PBS_QUEUE_DICT = [('bigdata', 24, 0, 24)]
    PBS_NJOB_LIMIT = 64

    GMX_MULTI = True
    GMX_MULTI_NJOB = 8
    GMX_MULTI_NOMP = 6  # set this if NGPU == 0


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
elif hostname.endswith('sjtu.edu.cn'):
    Config = PIConfig
elif hostname.startswith('ln') or hostname.startswith('cn'):
    Config = TH2Config
elif hostname == 'z-Mac.local':
    Config = MacConfig
else:
    raise Exception('msd-server will not work on this machine')
