import json, sys
from config import Config
sys.path.append(Config.MS_TOOLS_DIR)
from mstools.smiles.smiles import *


class ClassificationAtomType:
    def __init__(self, AtomicNum=[6]):
        self.AtomicNum = AtomicNum


class ClassificationSMILESSMARTS:
    def __init__(self, AtomicNum=[6]):
        self.AtomicNum = AtomicNum

    def classify(self, smiles):
        rdk_mol = Chem.MolFromSmiles(smiles)
        py_mol = pybel.readstring('smi', smiles)
        if not self.__AtomicNumCheck(rdk_mol, self.AtomicNum):
            return False
        else:
            return True

    @staticmethod
    def __AtomicNumCheck(rdk_mol, AtomicNum):
        for atom in rdk_mol.GetAtoms():
            if atom.GetAtomicNum() not in AtomicNum:
                return False
        return True

    @staticmethod
    def __has_ring(rdk_mol):
        return False if rdk_mol.GetRingInfo().NumRings() == 0 else True


def task_selection(task, select=True, rule='not select'):
    # basic_list: molecule must contain all atom types in basic list.
    # specific_list: molecule must contain at least one atom types in specific list.
    # addition_list：all the atom types in molecule must in (basic_list + specific_list + addition_list).
    # set addition_list = None to ignore this criterion
    # rule = 'ionic liquid'
    if rule == 'not select' or not select:
        return True
    elif rule == 'alkane':
        basic_list = ['h_1', 'c_4h3']
        specific_list = []
        addition_list = ['c_4h4', 'c_4h2', 'c_4h', 'c_4']
        reject_list = []
        ring = False
    elif rule == 'cyclic alkane':
        basic_list = ['h_1']
        specific_list = []
        addition_list = ['c_4h4', 'c_4h3', 'c_4h2', 'c_4h', 'c_4']
        reject_list = []
        ring = True
    elif rule == '3 ring':
        basic_list = ['h_1']
        specific_list = ['c_43', 'c_43h', 'c_43h2']
        addition_list = ['c_4h4', 'c_4h3', 'c_4h2', 'c_4h', 'c_4']
        reject_list = []
        ring = True
    elif rule == '4 ring':
        basic_list = ['h_1']
        specific_list = ['c_44', 'c_44h', 'c_44h2']
        addition_list = ['c_4h4', 'c_4h3', 'c_4h2', 'c_4h', 'c_4']
        reject_list = []
        ring = True
    elif rule == '5 ring':
        basic_list = ['h_1']
        specific_list = ['c_45', 'c_45h', 'c_45h2']
        addition_list = ['c_4h4', 'c_4h3', 'c_4h2', 'c_4h', 'c_4']
        reject_list = []
        ring = True
    elif rule == '5 ring alkene':
        basic_list = ['h_1']
        specific_list = ['c_35', 'c_35h']
        addition_list = ['c_4h4', 'c_4h3', 'c_4h2', 'c_4h', 'c_4', 'c_45', 'c_45h', 'c_45h2']
        reject_list = []
        ring = True
    elif rule == 'alkene':
        basic_list = ['h_1']
        specific_list = ['c_3', 'c_3h', 'c_3h2', 'c_3c3', 'c_3c3h']
        addition_list = ['c_4h4', 'c_4h3', 'c_4h2', 'c_4h', 'c_4']
        reject_list = []
        ring = False
    elif rule == 'cyclic alkene':
        basic_list = ['h_1']
        specific_list = ['c_3', 'c_3h', 'c_3h2', 'c_3c3', 'c_3c3h']
        addition_list = ['c_4h4', 'c_4h3', 'c_4h2', 'c_4h', 'c_4']
        reject_list = []
        ring = True
    elif rule == 'alkyne':
        basic_list = ['h_1']
        specific_list = ['c_2']
        addition_list = ['c_4h4', 'c_4h3', 'c_4h2', 'c_4h', 'c_4']
        reject_list = []
        ring = None
    elif rule == 'aromatics':
        basic_list = ['h_1']
        specific_list = ['c_3a', 'c_3ac']
        addition_list = ['c_4h4', 'c_4h3', 'c_4h2', 'c_4h', 'c_4']
        reject_list = []
        ring = None
    elif rule == 'CH':
        basic_list = ['h_1']
        specific_list = []
        addition_list = ['c_4h4', 'c_4h3', 'c_4h2', 'c_4h', 'c_4', 'c_3a', 'c_3ac', 'c_2', 'c_3', 'c_3h', 'c_3h2',
                         'c_3c3', 'c_3c3h', 'c_45', 'c_45h', 'c_45h2', 'c_44', 'c_44h', 'c_44h2', 'c_43', 'c_43h',
                         'c_43h2', 'c_35', 'c_35h']
        reject_list = []
        ring = None
    elif rule == 'alcohol CH':
        basic_list = ['h_1', 'h_1o', 'o_2h']
        specific_list = ['c_4oh3', 'c_4oh2', 'c_4oh', 'c_4o']
        addition_list = ['c_4h4', 'c_4h3', 'c_4h2', 'c_4h', 'c_4', 'c_3a', 'c_3ac', 'c_2', 'c_3', 'c_3h', 'c_3h2',
                         'c_3c3', 'c_3c3h', 'c_45', 'c_45h', 'c_45h2', 'c_44', 'c_44h', 'c_44h2', 'c_43', 'c_43h',
                         'c_43h2', 'c_35', 'c_35h']
        reject_list = []
        ring = None
    elif rule == 'ketone CH':
        basic_list = ['h_1', 'o_1']
        specific_list = ['c_3oh2', 'c_3oh', 'c_3o']
        addition_list = ['c_4h4', 'c_4h3', 'c_4h2', 'c_4h', 'c_4', 'c_3a', 'c_3ac', 'c_2', 'c_3', 'c_3h', 'c_3h2',
                         'c_3c3', 'c_3c3h', 'c_45', 'c_45h', 'c_45h2', 'c_44', 'c_44h', 'c_44h2', 'c_43', 'c_43h',
                         'c_43h2', 'c_35', 'c_35h']
        reject_list = []
        ring = None
    elif rule == 'alcohol':
        basic_list = ['h_1', 'h_1o', 'o_2h']
        specific_list = ['c_4oh3', 'c_4oh2', 'c_4oh', 'c_4o']
        addition_list = ['c_4h4', 'c_4h3', 'c_4h2', 'c_4h', 'c_4']
        reject_list = []
        ring = False
    elif rule == 'ketone':
        basic_list = ['h_1', 'o_1']
        specific_list = ['c_3oh2', 'c_3oh', 'c_3o']
        addition_list = ['c_4h4', 'c_4h3', 'c_4h2', 'c_4h', 'c_4']
        reject_list = []
        ring = False
    elif rule == 'ether':
        basic_list = ['h_1', 'o_2']
        specific_list = ['c_4oh3', 'c_4oh2', 'c_4oh', 'c_4o']
        addition_list = ['c_4h4', 'c_4h3', 'c_4h2', 'c_4h', 'c_4']
        reject_list = []
        ring = False
    elif rule == 'ammounia':
        basic_list = ['h_1']
        specific_list = ['n_3h2', 'n_3h', 'n_3', 'c_4nh3', 'c_4nh2', 'c_4nh', 'c_4n', 'h_1n']
        addition_list = ['c_4h4', 'c_4h3', 'c_4h2', 'c_4h', 'c_4']
        reject_list = []
        ring = False
    elif rule == 'acid':
        basic_list = ['h_1', 'o_1', 'o_2ch', 'h_1o']
        specific_list = ['c_3o2', 'c_3o2h']
        addition_list = ['c_4h4', 'c_4h3', 'c_4h2', 'c_4h', 'c_4']
        reject_list = []
        ring = False
    elif rule == 'imidazolium':
        basic_list = ['c_35an', 'c_35an2', 'c_35+da']
        specific_list = ['c_4np', 'c_4nph3']
        reject_list = []
        addition_list = None
        ring = None
    elif rule == 'ionic liquid':
        basic_list = []
        specific_list = []
        reject_list = ['br0-', 'i_0-']
        addition_list = None
        ring = None
    else:
        basic_list = []
        specific_list = []
        reject_list = []
        addition_list = None
        ring = None

    if task.atom_type is None:
        return False
    atom_type_list = json.loads(task.atom_type)
    for reject_atom in reject_list:
        if reject_atom in atom_type_list:
            return False
    if addition_list is not None:
        for atom_type in atom_type_list:
            if atom_type not in basic_list and atom_type not in specific_list and atom_type not in addition_list:
                return False
    for atom_type in basic_list:
        if atom_type not in atom_type_list:
            return False
    for atom_type in specific_list:
        if atom_type in atom_type_list:
            break
    else:
        if specific_list:
            return False
    smiles_list = json.loads(task.smiles_list)
    if ring is not None:
        ring_number = 0
        for smiles in smiles_list:
            ring_number += get_ring_number(smiles)
        if (ring and ring_number == 0) or (not ring and ring_number > 0):
            return False
    return True

'''
def IsCH(smiles):
    rdk_mol = Chem.MolFromSmiles(smiles)
    if rdk_mol is None:
        return False
    for atom in rdk_mol.GetAtoms():
        if atom.GetAtomicNum() not in [1, 6]:
            return False
    return True
'''