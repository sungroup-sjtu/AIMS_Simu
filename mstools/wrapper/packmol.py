import math
import os
import random
import subprocess
import sys
from subprocess import PIPE

from mstools.errors import PackmolError
from mstools.utils import count_atoms, greatest_common_divisor


class Packmol:
    '''
    wrappers for Packmol
    '''
    pass

    def __init__(self, packmol_bin):
        self.PACKMOL_BIN = packmol_bin

    def build_box(self, files: [str], numbers: [int], output: str,
                  natoms: int = None,
                  size: [float] = None, length: float = None,
                  tolerance: float = None,
                  seed: int = None) -> [int]:
        '''
        Build box directly from files

        :param files:
        :param numbers:
        :param output:
        :param size:
        :param length:
        :param tolerance:
        :param seed:
        :param save_input:
        :return:
        '''
        if len(files) == 0:
            raise PackmolError('No files provided')
        if len(files) != len(numbers):
            raise PackmolError('Invalid numbers')

        extensions = {filename.split('.')[-1].lower() for filename in files}
        if len(extensions) > 1:
            raise PackmolError('All file types should be the same')
        filetype = extensions.pop()

        if natoms != None:
            if natoms < 1:
                raise PackmolError('Invalid natoms')
            n_each_file = [count_atoms(filename) for filename in files]

            gcd_numbers = greatest_common_divisor(numbers)
            numbers = [int(i / gcd_numbers) for i in numbers]
            multiple = math.ceil(natoms / (sum([n_each_file[i] * number for i, number in enumerate(numbers)])))
            numbers = [multiple * i for i in numbers]

        if size != None:
            if len(size) != 3:
                raise PackmolError('Invalid box size')
            else:
                box = size
        elif length != None:
            box = [length, length, length]
        else:
            raise PackmolError('box size needed')

        tolerance = tolerance or 2.0

        seed = seed or random.randint(1e7, 1e8)

        inp = '''filetype {filetype}
tolerance {tolerance}
output {output}
seed {seed}
'''.format(filetype=filetype, tolerance=tolerance, output=output, seed=seed)

        for i, filename in enumerate(files):
            number = numbers[i]
            inp += '''
structure {filename}
  number {number}
  inside box 0 0 0 {box_size}
end structure
'''.format(filename=filename, number=number, box_size=' '.join(map(str, box)))

        with open('build.inp', 'w') as f:
            f.write(inp)

        # TODO subprocess PIPE not work on Mac, do not know why
        if sys.platform == 'darwin':
            os.system(self.PACKMOL_BIN + ' < build.inp > /dev/null')
        else:
            sp = subprocess.Popen([self.PACKMOL_BIN], stdin=PIPE, stdout=PIPE)
            sp.communicate(input=inp.encode())

        self.numbers = numbers
        self.length = box[0]
