from .jobmanager import JobManager


class Local(JobManager):
    def __init__(self):
        super().__init__()
        self.sh = '_job_local.sh'

    def generate_sh(self, workdir, commands, name=None):
        with open(self.sh, 'w') as f:
            f.write('#!/bin/sh\n\n'
                    'cd %(workdir)s\n'
                    % ({'workdir': workdir})
                    )
            for cmd in commands:
                f.write(cmd + '\n')

    def submit(self):
        print('Localhost is only for test')

    def get_info(self, name):
        raise Exception('Not supported on localhost')

    def kill_job(self, name):
        raise Exception('Not supported on localhost')
