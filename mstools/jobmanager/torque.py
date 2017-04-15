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
                        'workdir': workdir}
                       ))
            for cmd in commands:
                f.write(cmd + '\n')

    def submit(self):
        Popen(['qsub', self.sh]).communicate()
