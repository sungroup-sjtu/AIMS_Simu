import math
import os

from typing import List

from .errors import OpenBabelError


def count_atoms(filename):
    '''
    count atom numbers in PDB file
    '''
    if not os.path.exists(filename):
        raise Exception('file not exist')
    filetype = filename.split('.')[-1].lower()
    if filetype == 'pdb':
        with open(filename) as f:
            strf = f.read()
            n = strf.count('ATOM')
            if n == 0:
                n = strf.count('HETATM')
            return n
    else:
        return 1


def greatest_common_divisor(numbers):
    '''
    calculte the greatest common divisor
    '''
    minimal = min(numbers)
    for i in range(minimal, 1, -1):
        flag = True
        for number in numbers:
            if number % i != 0:
                flag = False
        if flag:
            return i
    return 1


def random_string(length=8):
    import random, string
    return ''.join(random.sample(string.ascii_letters, length))


def cd_or_create_and_cd(dir):
    if not os.path.exists(dir):
        try:
            os.mkdir(dir)
        except:
            raise Exception('Cannot create directory: %s' % dir)

    try:
        os.chdir(dir)
    except:
        raise Exception('Cannot read directory: %s' % dir)


def get_T_list_from_range(t_min, t_max, interval=20) -> List[int]:
    T_list = []

    t_min_new = int(t_min / interval)
    t_max_new = math.ceil(t_max / interval)
    for i in range(t_min_new, t_max_new + 1):
        T_list.append(i * interval)
    return T_list


def get_P_list_from_range(p_min, p_max, multiple=(5,)) -> List[int]:
    P_list = []

    multiple = list(multiple)
    if 1 not in multiple:
        multiple.append(1)
    multiple.sort()

    magnitude_min = int(math.log10(p_min))
    magnitude_max = math.ceil(math.log10(p_max))

    for i in range(magnitude_min, magnitude_max):
        for m in multiple:
            P_list.append(10 ** i * m)
    P_list.append(10 ** magnitude_max)
    return P_list


def create_pdb_from_smiles(smiles: str, pdb_out: str):
    try:
        import pybel
        py_mol = pybel.readstring('smi', smiles)
        py_mol.addh()
        py_mol.make3D()
        py_mol.write('pdb', pdb_out, overwrite=True)
    except:
        raise OpenBabelError('Cannot create PDB from SMILES')
    else:
        return py_mol
