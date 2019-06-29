import os
import sys

CWD = os.path.dirname(os.path.abspath(__file__))

class Config:
    SQLALCHEMY_BINDS = {
        'cv'  : 'sqlite:///%s?check_same_thread=False' % os.path.join(CWD, 'database/cv.sqlite'),
        'yaws': 'sqlite:///%s?check_same_thread=False' % os.path.join(CWD, 'database/yaws.sqlite'),
        'nist': 'sqlite:///%s?check_same_thread=False' % os.path.join(CWD, 'database/nist.sqlite'),
        'ilthermo': 'sqlite:///%s?check_same_thread=False' % os.path.join(CWD, 'database/ilthermo.sqlite'),
    }
    SQLALCHEMY_COMMIT_ON_TEARDOWN = True
    SQLALCHEMY_TRACK_MODIFICATIONS = True
    
    
    MS_TOOLS_DIR = os.path.join(CWD, '..', 'AIMS_Tools')
    WORK_DIR = os.path.join(CWD, 'SimulationData')
    RUN_DIR = os.path.join(CWD, 'run')

    PACKMOL_BIN = '/share/apps/tools/packmol'
    DFF_ROOT = '/share/workspace/xiangyan/src/DFF/Developing'
    # DFF_TABLE = 'IL'
    # DFF_TABLE = 'MGI'

    CHARGE_SCALE = 0.8 if DFF_TABLE == 'IL' else 1.0
    PPF_FILE = os.path.join(RUN_DIR, 'ff.ppf')
    PPF = PPF_FILE if os.path.exists(PPF_FILE) else None

    # simulation details setting
    NATOMS = 3000 # least number of atoms build in simulation box.
    NMOLS = 120 # least number of molecules build in simulation box.
    LJ96 = False # using LJ 9-6 non-bonded potential
    DIFF_GK = False # using green-kubo method to calculate the diffusion constant. (Expensive, not suggest)
    DEBUG = False # if true: do not delete the trajectory file in analyze process.

    EXTEND_CYCLE_LIMIT = 20


    @classmethod
    def init_app(cls, app):
        sys.path.append(Config.MS_TOOLS_DIR)
        from mstools.jobmanager import Slurm, Torque, RemoteSlurm

        _pbs_dict = {'slurm': Slurm, 'remote_slurm': RemoteSlurm, 'torque': Torque}
        # run
        PBS = _pbs_dict[cls.PBS_MANAGER]
        jobmanager = PBS(*cls.PBS_ARGS, **cls.PBS_KWARGS)
        if hasattr(cls, 'PBS_SUBMIT_CMD'):
            jobmanager.submit_cmd = cls.PBS_SUBMIT_CMD
        if hasattr(cls, 'PBS_TIME_LIMIT'):
            jobmanager.time = cls.PBS_TIME_LIMIT
        # extend
        PBS = _pbs_dict[cls.EXTEND_PBS_MANAGER]
        jm_extend = PBS(*cls.EXTEND_PBS_ARGS, **cls.EXTEND_PBS_KWARGS)
        if hasattr(cls, 'EXTEND_PBS_SUBMIT_CMD'):
            jm_extend.submit_cmd = cls.EXTEND_PBS_SUBMIT_CMD
        if hasattr(cls, 'EXTEND_PBS_TIME_LIMIT'):
            jm_extend.time = cls.EXTEND_PBS_TIME_LIMIT

        if jm_extend.is_remote:
            raise Exception('Remote jobmanager is not compatible with extend')

        # bugfix
        PBS = _pbs_dict[cls.BUGFIX_PBS_MANAGER]
        jm_bugfix = PBS(*cls.BUGFIX_PBS_ARGS, **cls.BUGFIX_PBS_KWARGS)
        if hasattr(cls, 'BUGFIX_PBS_SUBMIT_CMD'):
            jm_bugfix.submit_cmd = cls.BUGFIX_PBS_SUBMIT_CMD
        if hasattr(cls, 'BUGFIX_PBS_TIME_LIMIT'):
            jm_bugfix.time = cls.BUGFIX_PBS_TIME_LIMIT

        if jm_bugfix.is_remote:
            raise Exception('Remote jobmanager is not compatible with bugfix')


        app.jobmanager = jobmanager
        app.jm_extend = jm_extend
        app.jm_bugfix = jm_bugfix


class SunRunConfig:
    PBS_NJOB_LIMIT = 100
    PBS_MANAGER = 'slurm'
    # Use GPU 
    PBS_ARGS = ('gtx', 32, 2, 16)  # partition, cpu(hyperthreading), gpu, cpu_request
    GMX_MDRUN= 'gmx_gpu mdrun'
    GMX_MULTI = True
    GMX_MULTI_NJOB = 8 # Use -multidir function of GROMACS. For Npt simulation, set it to 8. For NvtSlab simulation, 4 is better
    GMX_MULTI_NOMP = None  # Set the OpenMP threads. When set to None, use only one node and the best number of threads is automatically determined

    # Use CPU
    # PBS_ARGS = ('cpu', 8, 0, 8)  # partition, cpu(hyperthreading), gpu, cpu_request
    # GMX_MDRUN= 'gmx_gpu mdrun'
    # GMX_MULTI = False
    # GMX_MULTI_NJOB = 8 # Use -multidir function of GROMACS. For Npt simulation, set it to 8. For NvtSlab simulation, 4 is better
    # GMX_MULTI_NOMP = None  # Set the OpenMP threads. When set to None, use only one node and the best number of threads is automatically determined

    # Use Fast
    # PBS_ARGS = ('fast', 12, 0, 12)  # partition, cpu(hyperthreading), gpu, cpu_request
    # GMX_MDRUN= 'gmx_fast mdrun'
    # GMX_MULTI = True
    # GMX_MULTI_NJOB = 4 # Use -multidir function of GROMACS. For Npt simulation, set it to 8. For NvtSlab simulation, 4 is better
    # GMX_MULTI_NOMP = None  # Set the OpenMP threads. When set to None, use only one node and the best number of threads is automatically determined

    PBS_KWARGS = {'env_cmd': 'module purge; module load icc gromacs/2018.6'}
    PBS_TIME_LIMIT = 240  # hour
    GMX_BIN = '/share/apps/gromacs/2018.6/bin/gmx_serial'
    

class SunExtendConfig:
    EXTEND_PBS_NJOB_LIMIT = 200
    EXTEND_PBS_MANAGER = 'slurm'
    # Use GPU
    EXTEND_PBS_ARGS = ('gtx', 32, 2, 16)
    EXTEND_GMX_MDRUN = 'gmx_gpu mdrun'
    EXTEND_GMX_MULTI = True
    EXTEND_GMX_MULTI_NJOB = 8

    # Use CPU
    # EXTEND_PBS_ARGS = ('cpu', 8, 0, 8)
    # EXTEND_GMX_MDRUN = 'gmx_gpu mdrun'
    # EXTEND_GMX_MULTI = False #

    # Use FAST
    # EXTEND_PBS_ARGS = ('fast', 12, 0, 12)
    # EXTEND_GMX_MDRUN = 'gmx_fast mdrun'
    # EXTEND_GMX_MULTI = False #

    EXTEND_PBS_KWARGS = {'env_cmd': 'module purge; module load icc gromacs/2018.6'}
    EXTEND_PBS_TIME_LIMIT = 240  # hour

    EXTEND_GMX_BIN = '/share/apps/gromacs/2018.6/bin/gmx_serial'


class SunBugFixConfig:
    BUGFIX_PBS_NJOB_LIMIT = 200
    BUGFIX_PBS_MANAGER = 'slurm'
    # Use GPU
    # BUGFIX_PBS_ARGS = ('gtx', 32, 2, 16)
    # BUGFIX_GMX_MDRUN = 'gmx_gpu mdrun'
    # BUGFIX_GMX_MULTI = True
    # BUGFIX_GMX_MULTI_NJOB = 8

    # Use CPU
    # BUGFIX_PBS_ARGS = ('cpu', 8, 0, 8)
    # BUGFIX_GMX_MDRUN = 'gmx_gpu mdrun'
    # BUGFIX_GMX_MULTI = False #

    # Use FAST
    BUGFIX_PBS_ARGS = ('fast', 12, 0, 12)
    BUGFIX_GMX_MDRUN = 'gmx_fast mdrun'
    BUGFIX_GMX_MULTI = False #

    BUGFIX_PBS_KWARGS = {'env_cmd': 'module purge; module load icc gromacs/2018.6'}
    BUGFIX_PBS_TIME_LIMIT = 240  # hour

    BUGFIX_GMX_BIN = '/share/apps/gromacs/2018.6/bin/gmx_serial'


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


class NptConfig(Config, SunRunConfig, SunExtendConfig, SunBugFixConfig):
    DB = 'database/msd.npt.db'
    SQLALCHEMY_DATABASE_URI = 'sqlite:///%s?check_same_thread=False' % os.path.join(CWD, DB)
    LOG = os.path.join(CWD, '_log.npt.txt')
    ALLOWED_PROCEDURES = ['npt', 'ppm', 'npt-v-rescale', 'npt-2', 'npt-3']
    REPEAT_NUMBER = 1

class NptMultiConfig(Config, SunRunConfig, SunExtendConfig, SunBugFixConfig):
    DB = 'database/msd.npt.db'
    SQLALCHEMY_DATABASE_URI = 'sqlite:///%s?check_same_thread=False' % os.path.join(CWD, DB)
    LOG = os.path.join(CWD, '_log.npt_multi.txt')
    ALLOWED_PROCEDURES = ['npt-multi']
    REPEAT_NUMBER = 80

class NvtMultiConfig(Config, SunRunConfig, SunExtendConfig, SunBugFixConfig):
    DB = 'database/msd.npt.db'
    SQLALCHEMY_DATABASE_URI = 'sqlite:///%s?check_same_thread=False' % os.path.join(CWD, DB)
    LOG = os.path.join(CWD, '_log.nvt_multi.txt')
    ALLOWED_PROCEDURES = ['nvt-multi', 'nvt-multi-2', 'nvt-multi-3']
    REPEAT_NUMBER = 80

class NvtSlabConfig(Config, SunRunConfig, SunExtendConfig, SunBugFixConfig):
    DB = 'database/msd.nvt-slab.db'
    SQLALCHEMY_DATABASE_URI = 'sqlite:///%s?check_same_thread=False' % os.path.join(CWD, DB)
    LOG = os.path.join(CWD, '_log.nvt-slab.txt')
    ALLOWED_PROCEDURES = ['nvt-slab']
    REPEAT_NUMBER = 1


configs = {
    'npt': NptConfig,
    'npt-2': NptConfig,
    'npt-3': NptConfig,
    'npt-v-rescale': NptConfig,
    'npt-multi': NptMultiConfig,
    'nvt-slab': NvtSlabConfig,
    'ppm': NptConfig,
    'nvt-multi': NvtMultiConfig,
    'nvt-multi-2': NvtMultiConfig,
    'nvt-multi-3': NvtMultiConfig,
}
