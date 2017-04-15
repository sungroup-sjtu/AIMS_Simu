import math

from ..utils import create_pdb_from_smiles
from ..wrapper import Packmol, DFF
from .procedure import Procedure
from ..jobmanager import Local


class Simulation():
    def __init__(self, packmol_bin=None, dff_root=None, procedure=None, jobmanager=None):
        self.packmol = Packmol(packmol_bin=packmol_bin)
        self.dff = DFF(dff_root=dff_root)
        self.procedure = procedure
        self.jobmanager = jobmanager

        self.n_mol: int
        self.msd: str = 'init.msd'

    def build(self, smiles, n_atoms):
        pass

    def prepare(self, model_dir):
        pass

    def run(self):
        self.jobmanager.submit()

    def analyze(self):
        pass

    def build_dff_box_from_smiles(self, smiles, n_atoms):
        py_mol = create_pdb_from_smiles(smiles, 'mol.pdb')
        n_atom_per_mol = len(py_mol.atoms)
        self.n_mol = math.ceil(n_atoms / n_atom_per_mol)  # A total of 3000 atoms
        # length = (10 / 6.022 * py_mol.molwt * n_mol / density) ** (1 / 3) # mass density
        length = (9.96 * self.n_mol * n_atom_per_mol) ** (1 / 3)  # number density according to water

        print('Build coordinates using Packmol: %i molecules ...' % self.n_mol)
        if self.procedure == Procedure.NPT_BINARY_SLAB:
            self.packmol.build_box(['mol.pdb'], [self.n_mol], 'init.pdb', length=length - 2, slab=True, silent=True)
        else:
            self.packmol.build_box(['mol.pdb'], [self.n_mol], 'init.pdb', length=length - 2, silent=True)

        print('Create box using DFF ...')
        if self.procedure == Procedure.NVT_SLAB:
            self.dff.build_bulk_after_packmol('mol.pdb', self.n_mol, self.msd, pdb_corr='init.pdb',
                                              size=[length, length, length * 3])
        else:
            self.dff.build_bulk_after_packmol('mol.pdb', self.n_mol, self.msd, pdb_corr='init.pdb', length=length)
