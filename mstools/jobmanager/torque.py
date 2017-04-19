import subprocess
from subprocess import Popen

from .jobmanager import JobManager


class Torque(JobManager):
    def __init__(self, queue, nprocs):
        super().__init__(queue=queue, nprocs=nprocs)
        self.sh = 'job_torque.sh'
        self.out = 'job_torque.out'
        self.err = 'job_torque.err'

    def generate_sh(self, workdir, commands, name):
        with open(self.sh, 'w') as f:
            f.write('#!/bin/sh\n'
                    '#PBS -N %(name)s\n'
                    '#PBS -o %(out)s\n'
                    '#PBS -e %(err)s\n'
                    '#PBS -q %(queue)s\n'
                    '#PBS -l nodes=1:ppn=%(nprocs)s\n\n'
                    'cd %(workdir)s\n'
                    % ({'name': name,
                        'out': self.out,
                        'err': self.err,
                        'queue': self.queue,
                        'nprocs': self.nprocs,
                        'workdir': workdir
                        })
                    )
            for cmd in commands:
                f.write(cmd + '\n')

    def submit(self):
        Popen(['qsub', self.sh]).communicate()

    def get_info_from_id(self, id) -> bool:
        try:
            output = subprocess.check_output(['qstat', '-f', str(id)])
        except:
            return False

        for line in output.decode().splitlines():
            if line.strip().startswith('job_state'):
                state = line.split()[-1]
                if state in ['R', 'Q']:
                    return True
                else:
                    return False
        return False

    def get_id_from_name(self, name: str) -> int:
        try:
            output = subprocess.check_output(['qstat'])
        except:
            raise

        for line in output.decode().splitlines():
            if line.find(name) != -1:
                return int(line.strip().split()[0])

        return None

    def get_info_from_name(self, name) -> bool:
        id = self.get_id_from_name(name)
        if id == None:
            return False
        return self.get_info_from_id(id)
