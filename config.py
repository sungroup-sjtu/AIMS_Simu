import os, socket


class Config:
    DB_NAME = 'msdserver.sqlite'
    DB_PATH = os.path.dirname(os.path.abspath(__file__)) + os.sep + DB_NAME

    SQLALCHEMY_DATABASE_URI = 'sqlite:///%s' % DB_PATH
    SQLALCHEMY_COMMIT_ON_TEARDOWN = True
    SQLALCHEMY_TRACK_MODIFICATIONS = True

    if socket.gethostname() == 'cluster.hpc.org':
        USER = 'msdserver'
        WORK_DIR = '/share/workspace/msdserver/MSDataServer/'
        MS_TOOLS_PATH = '/home/gongzheng/GitHub/ms-tools'
        DFF_ROOT = '/share/apps/DFF/7.2'
        PACKMOL_BIN = '/share/apps/tools/packmol'
        LAMMPS_BIN = '/share/apps/lammps/lmp-stable'
        GMX_BIN = '/share/apps/gromacs/2016.2/bin/gmx'
    elif socket.gethostname() == 'z-Mac.local':
        USER = 'zheng'
        WORK_DIR = '/tmp/MSDataServer/'
        MS_TOOLS_PATH = '/Users/zheng/Projects/ms-tools'
        DFF_ROOT = '/Users/zheng/Projects/DFF7.2'
        PACKMOL_BIN = '/Users/zheng/Projects/DFF7.2/bin32m/Packmol/packmol.exe'
        LAMMPS_BIN = '/usr/local/bin/lmp_serial'
        GMX_BIN = '/usr/local/bin/gmx'
    else:
        raise Exception('MSData Server will not work on this machine')
