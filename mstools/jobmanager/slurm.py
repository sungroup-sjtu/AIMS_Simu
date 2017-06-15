import subprocess
from subprocess import Popen

from .jobmanager import JobManager


class Slurm(JobManager):
    def __init__(self, queue, nprocs):
        super().__init__(queue=queue, nprocs=nprocs)
        self.sh = '_job_slurm.sh'
        self.out = '_job_slurm.out'
        self.err = '_job_slurm.err'

    def generate_sh(self, workdir, commands, name):
        with open(self.sh, 'w') as f:
            f.write('#!/bin/bash\n'
                    '#SBATCH --job-name=%(name)s\n'
                    '#SBATCH --partition=%(queue)s\n'
                    '#SBATCH --output=%(out)s\n'
                    '#SBATCH --error=%(err)s\n'
                    '#SBATCH -n %(nprocs)s\n'
                    '#SBATCH --tasks-per-node=%(tasks)s\n\n'
                    'source /usr/share/Modules/init/bash\n'
                    'unset MODULEPATH\n'
                    'module use /lustre/usr/modulefiles/pi\n'
                    'module purge\n'
                    'module load icc/16.0 impi/5.1 mkl/11.3\n'
                    'export I_MPI_PMI_LIBRARY=/usr/lib64/libpmi.so\n'
                    'export I_MPI_FABRICS=shm:dapl\n\n'
                    'source $HOME/software/gromacs/5.1.4-msdserver/bin/GMXRC.bash\n\n'
                    'cd %(workdir)s\n\n'
                    % ({'name': name,
                        'out': self.out,
                        'err': self.err,
                        'queue': self.queue,
                        'nprocs': self.nprocs,
                        'tasks': min(self.nprocs, 16),
                        'workdir': workdir
                        })
                    )
            for cmd in commands:
                f.write(cmd + '\n')

    def submit(self):
        Popen(['sbatch', self.sh]).communicate()

    def get_info_from_id(self, id) -> bool:
        #try:
            #output = subprocess.check_output(['qstat', '-f', str(id)])
        #except:
            #return False

        #for line in output.decode().splitlines():
            #if line.strip().startswith('job_state'):
                #state = line.split()[-1]
                #if state in ['R', 'Q']:
                    #return True
                #else:
                    #return False
        return False

    def get_id_from_name(self, name: str) -> int:
        #try:
            #output = subprocess.check_output(['qstat'])
        #except:
            #raise

        #for line in output.decode().splitlines():
            #if line.find(name) != -1:
                #return int(line.strip().split()[0])

        return None

    def get_info(self, name) -> bool:
        id = self.get_id_from_name(name)
        if id == None:
            return False
        return self.get_info_from_id(id)

    def kill_job(self, name) -> bool:
        #id = self.get_id_from_name(name)
        #if id == None:
            #return False
        #try:
            #subprocess.check_call(['qdel', str(id)])
        #except:
            #raise Exception('Cannot kill job: %s' % name)

        return True
