from subprocess import Popen

from .jobmanager import JobManager


class Local(JobManager):
    def __init__(self):
        super().__init__()
        self.sh = 'job_local.sh'

    def generate_sh(self, workdir, commands, name=None):
        with open(self.sh, 'w') as f:
            f.write('#!/bin/sh\n\n'
                    'cd %(workdir)s\n'
                    % (workdir))
            for cmd in commands:
                f.write(cmd + '\n')

    def submit(self):
        Popen(['sh', self.sh]).communicate()
