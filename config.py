import os
import sys

CWD = os.path.dirname(os.path.abspath(__file__))


class Config:
    SQLALCHEMY_BINDS = {
        'cv'  : 'sqlite:///%s?check_same_thread=False' % os.path.join(CWD, 'database/cv.sqlite'),
        'yaws': 'sqlite:///%s?check_same_thread=False' % os.path.join(CWD, 'database/yaws.sqlite'),
        'nist': 'sqlite:///%s?check_same_thread=False' % os.path.join(CWD, 'database/nist.sqlite'),
    }
    SQLALCHEMY_COMMIT_ON_TEARDOWN = True
    SQLALCHEMY_TRACK_MODIFICATIONS = True

    MS_TOOLS_DIR = os.path.join(CWD, '../ms-tools')

    WORK_DIR = '/share/md1400/_MSDServer/'
    PACKMOL_BIN = '/share/apps/tools/packmol'
    DFF_ROOT = '/home/gongzheng/apps/DFF/Developing'
    DFF_TABLE = 'MGI'

    EXTEND_CYCLE_LIMIT = 8

    @classmethod
    def init_app(cls, app):
        sys.path.append(Config.MS_TOOLS_DIR)
        from mstools.jobmanager import Slurm, Torque, RemoteSlurm

        _pbs_dict = {'slurm': Slurm, 'remote_slurm': RemoteSlurm, 'torque': Torque}

        PBS = _pbs_dict[cls.PBS_MANAGER]
        jobmanager = PBS(*cls.PBS_ARGS, **cls.PBS_KWARGS)
        if hasattr(cls, 'PBS_SUBMIT_CMD'):
            jobmanager.submit_cmd = cls.PBS_SUBMIT_CMD
        if hasattr(cls, 'PBS_TIME_LIMIT'):
            jobmanager.time = cls.PBS_TIME_LIMIT

        PBS = _pbs_dict[cls.EXTEND_PBS_MANAGER]
        jm_extend = PBS(*cls.EXTEND_PBS_ARGS, **cls.EXTEND_PBS_KWARGS)
        if hasattr(cls, 'EXTEND_PBS_SUBMIT_CMD'):
            jm_extend.submit_cmd = cls.EXTEND_PBS_SUBMIT_CMD
        if hasattr(cls, 'EXTEND_PBS_TIME_LIMIT'):
            jm_extend.time = cls.EXTEND_PBS_TIME_LIMIT

        if jm_extend.is_remote:
            raise Exception('Remote jobmanager is not compatible with extend')

        app.jobmanager = jobmanager
        app.jm_extend = jm_extend


class SunRunConfig:
    PBS_NJOB_LIMIT = 24
    PBS_MANAGER = 'slurm'
    # PBS_ARGS = ('cpu', 64, 0, 32)  # partition, cpu(hyperthreading), gpu, cpu_request
    PBS_ARGS = ('gtx', 32, 2, 16)  # partition, cpu(hyperthreading), gpu, cpu_request
    # PBS_ARGS = ('fast', 24, 0, 12)  # partition, cpu(hyperthreading), gpu, cpu_request
    PBS_KWARGS = {'env_cmd': 'module purge; module load gcc openmpi gromacs/2016.5'}
    PBS_TIME_LIMIT = 10  # hour

    GMX_BIN = 'gmx_gpu'
    # GMX_BIN = 'gmx_fast'
    GMX_MDRUN = None
    GMX_MULTI = True
    GMX_MULTI_NJOB = 4
    GMX_MULTI_NOMP = None  # Use only one node. Automatically determine the best number of threads.


class SunExtendConfig:
    EXTEND_PBS_NJOB_LIMIT = 24
    EXTEND_PBS_MANAGER = 'slurm'
    EXTEND_PBS_ARGS = ('cpu', 64, 0, 32)
    # EXTEND_PBS_ARGS = ('gtx', 32, 2, 16)
    # EXTEND_PBS_ARGS = ('fast', 24, 0, 12)
    EXTEND_PBS_KWARGS = {'env_cmd': 'module purge; module load gcc openmpi gromacs/2016.5'}
    EXTEND_PBS_TIME_LIMIT = 10  # hour

    EXTEND_GMX_BIN = 'gmx_gpu'
    # EXTEND_GMX_BIN = 'gmx_fast'
    EXTEND_GMX_MDRUN = None
    EXTEND_GMX_MULTI = False
    EXTEND_GMX_MULTI_NJOB = 2


class PiRunConfig:
    _env_cmd = '''
source /usr/share/Modules/init/bash
module purge; unset MODULEPATH; module use /lustre/usr/modulefiles; module load icc/16.0 mkl/2016 impi/2016
export I_MPI_PMI_LIBRARY=/usr/lib64/libpmi.so; export I_MPI_FABRICS=shm:dapl
'''

    PBS_NJOB_LIMIT = 180
    PBS_MANAGER = 'remote_slurm'
    PBS_ARGS = ('cpu', 16, 0, 16,)  # partition, cpu(hyperthreading), gpu, cpu_request
    PBS_KWARGS = {
        'host'      : '202.120.58.230',
        'username'  : 'nishsun-1',
        'remote_dir': '/lustre/home/acct-nishsun/nishsun-1/workspace/_REMOTE_SLURM_/',
        'env_cmd'   : _env_cmd
    }
    # PBS_SUBMIT_CMD = 'sbatch --reservation=cpu_nishsun'

    GMX_BIN = 'gmx_serial'
    GMX_MDRUN = 'gmx mdrun'
    GMX_MULTI = True
    GMX_MULTI_NJOB = 2
    GMX_MULTI_NOMP = 8


class NptConfig(Config, SunRunConfig, SunExtendConfig):
    DB = 'database/msd.npt.db'
    SQLALCHEMY_DATABASE_URI = 'sqlite:///%s?check_same_thread=False' % os.path.join(CWD, DB)
    LOG = os.path.join(CWD, '_log.npt.txt')
    ALLOWED_PROCEDURES = ['npt']


class NvtSlabConfig(Config, SunRunConfig, SunExtendConfig):
    DB = 'database/msd.nvt-slab.db'
    SQLALCHEMY_DATABASE_URI = 'sqlite:///%s?check_same_thread=False' % os.path.join(CWD, DB)
    LOG = os.path.join(CWD, '_log.nvt-slab.txt')
    ALLOWED_PROCEDURES = ['nvt-slab']


configs = {
    'npt'     : NptConfig,
    'nvt-slab': NvtSlabConfig,
}
