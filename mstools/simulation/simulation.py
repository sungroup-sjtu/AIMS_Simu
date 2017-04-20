import math

from ..utils import create_mol_from_smiles
from ..wrapper import Packmol, DFF


class Simulation():
    def __init__(self, packmol_bin=None, dff_root=None, jobmanager=None):
        self.packmol = Packmol(packmol_bin=packmol_bin)
        self.dff = DFF(dff_root=dff_root)
        self.jobmanager = jobmanager
        self.procedure = None

        self.n_mol_list: [int]
        self.msd: str = 'init.msd'

    def set_procedure(self, procedure):
        self.procedure = procedure

    def build(self):
        pass

    def prepare(self, model_dir):
        pass

    def run(self):
        self.jobmanager.submit()

    def check_finished(self):
        pass

    def analyze(self) -> (bool, {str: float}):
        pass

    def set_system(self, smiles_list: [str], n_atoms: int):
        self.pdb_list = []
        self.mol2_list = []
        n_components = len(smiles_list)
        n_atom_per_mol_list = []
        for i, smiles in enumerate(smiles_list):
            pdb = 'mol-%i.pdb' % (i + 1)
            mol2 = 'mol-%i.mol2' % (i + 1)
            py_mol = create_mol_from_smiles(smiles, pdb_out=pdb, mol2_out=mol2)
            self.pdb_list.append(pdb)
            self.mol2_list.append(mol2)
            n_atom_per_mol_list.append(len(py_mol.atoms))

        self.n_mol_list = [math.ceil(n_atoms / n_components / n_atom) for n_atom in n_atom_per_mol_list]
        n_atoms = sum([n_atom_per_mol_list[i] * self.n_mol_list[i] for i in range(n_components)])
        # length = (10 / 6.022 * py_mol.molwt * n_mol / density) ** (1 / 3) # mass density
        density_relative_to_water = 0.8
        self.length = (9.96 * n_atoms / density_relative_to_water) ** (1 / 3)  # number density according to water
