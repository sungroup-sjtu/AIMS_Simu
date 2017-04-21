import math
import os

from .errors import OpenBabelError


def greatest_common_divisor(numbers):
    '''
    calculate the greatest common divisor
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


def get_T_list_from_range(t_min: int, t_max: int, interval: int = None, number: int = None) -> [int]:
    t_min = max(1, t_min)
    t_max = max(t_min, t_max)
    if t_min == t_max:
        return [t_min]

    if number != None:
        if number <= 1:
            interval = t_max - t_min
        else:
            interval = math.ceil((t_max - t_min) / (number - 1))
    elif interval == None:
        interval = 20

    interval = max(5, interval)

    if interval > t_max - t_min:
        return [t_min, t_max]

    T_list = []
    i = 0
    while True:
        t = t_min + i * interval
        T_list.append(t)
        if t >= t_max:
            break
        i += 1

    return T_list


def get_P_list_from_range(p_min, p_max, multiple=(5,), n_point: int = None) -> [int]:
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


def create_mol_from_smiles(smiles: str, pdb_out: str = None, mol2_out: str = None):
    try:
        import pybel
        py_mol = pybel.readstring('smi', smiles)
        py_mol.addh()
        py_mol.make3D()
        if pdb_out != None:
            py_mol.write('pdb', pdb_out, overwrite=True)
        if mol2_out != None:
            py_mol.write('mol2', mol2_out, overwrite=True)
    except:
        raise OpenBabelError('Cannot create molecule from SMILES')
    else:
        return py_mol


def estimate_density_from_formula(f) -> float:
    from .formula import Formula
    # unit: g/mL
    return Formula.read(f).estimate_density()


def check_converged(data_series: [float]) -> bool:
    return True
