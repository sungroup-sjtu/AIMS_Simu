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
    '''
    This config determines the PBS and GROMACS information for running simulation. Should be very careful
    GMX_BIN is running on the main node. The environment variables are configured by .bashrc or .zshrc or whatever
    GMX_MDRUN is running on the compute node. The environment variables are configured by PBS_KWARGS['env_cmd']
    Make sure that the version of GMX_BIN and GMX_MDRUN are the same
    '''
    PBS_NJOB_LIMIT = 48
    PBS_MANAGER = 'slurm'
    # PBS_ARGS = ('cpu', 64, 0, 32)  # partition, cpu(hyperthreading), gpu, cpu_request
    PBS_ARGS = ('gtx', 32, 2, 16)  # partition, cpu(hyperthreading), gpu, cpu_request
    # PBS_ARGS = ('fast', 24, 0, 12)  # partition, cpu(hyperthreading), gpu, cpu_request
    PBS_KWARGS = {'env_cmd': 'module purge; module load icc gromacs/2016.6'}
    PBS_TIME_LIMIT = 1  # hour

    GMX_BIN = '/share/apps/gromacs/2016.6/bin/gmx_serial'
    GMX_MDRUN = 'gmx_gpu mdrun'
    # GMX_MDRUN= 'gmx_fast mdrun'
    GMX_MULTI = True
    GMX_MULTI_NJOB = 8  # Use -multidir function of GROMACS. For Npt simulation, set it to 8. For NvtSlab simulation, 4 is better
    GMX_MULTI_NOMP = None  # Set the OpenMP threads. When set to None, use only one node and the best number of threads is automatically determined


class SunExtendConfig:
    '''
    This config determines the PBS and GROMACS information for extending simulation
    '''
    EXTEND_PBS_NJOB_LIMIT = 48
    EXTEND_PBS_MANAGER = 'slurm'
    # EXTEND_PBS_ARGS = ('cpu', 64, 0, 32)
    # EXTEND_PBS_ARGS = ('gtx', 32, 2, 16)
    EXTEND_PBS_ARGS = ('fast', 24, 0, 12)
    EXTEND_PBS_KWARGS = {'env_cmd': 'module purge; module load icc gromacs/2016.6'}
    EXTEND_PBS_TIME_LIMIT = 1  # hour

    EXTEND_GMX_BIN = '/share/apps/gromacs/2016.6/bin/gmx_serial'
    # EXTEND_GMX_MDRUN = 'gmx_gpu mdrun'
    EXTEND_GMX_MDRUN = 'gmx_fast mdrun'
    EXTEND_GMX_MULTI = False  # Do not user -multidir function for extend in case of weird bugs
    EXTEND_GMX_MULTI_NJOB = 2


class PiRunConfig:
    _env_cmd = '''
source /usr/share/Modules/init/bash
module purge; unset MODULEPATH; module use /lustre/usr/modulefiles; module load icc/16.0 mkl/2016 impi/2016
export I_MPI_PMI_LIBRARY=/usr/lib64/libpmi.so; export I_MPI_FABRICS=shm:dapl
'''

    PBS_NJOB_LIMIT = 180
    PBS_MANAGER = 'remote_slurm' # Use slurm on remote hpc. Do not use it unless you know exactly the mechanism
    PBS_ARGS = ('cpu', 16, 0, 16,)  # partition, cpu(hyperthreading), gpu, cpu_request
    PBS_KWARGS = {
        'host'      : os.getenv('PI_HOST'),
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
