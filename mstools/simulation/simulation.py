import math

from .procedure import Procedure
from ..utils import create_pdb_from_smiles
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

    def build(self, smiles_list: [str], n_atoms: int):
        pass

    def prepare(self, model_dir):
        pass

    def run(self):
        self.jobmanager.submit()

    def check_finished(self):
        pass

    def analyze(self):
        pass

    def build_dff_box_from_smiles(self, smiles_list: [str], n_atoms: int):
        n_components = len(smiles_list)
        pdb_list = []
        n_atom_per_mol_list = []
        for i, smiles in enumerate(smiles_list):
            pdb = 'mol-%i.pdb' % (i + 1)
            py_mol = create_pdb_from_smiles(smiles, pdb)
            pdb_list.append(pdb)
            n_atom_per_mol_list.append(len(py_mol.atoms))

        self.n_mol_list = [math.ceil(n_atoms / n_components / n_atom) for n_atom in n_atom_per_mol_list]
        n_atoms = sum([n_atom_per_mol_list[i] * self.n_mol_list[i] for i in range(n_components)])
        # length = (10 / 6.022 * py_mol.molwt * n_mol / density) ** (1 / 3) # mass density
        density_relative_to_water = 0.8
        length = (9.96 * n_atoms / density_relative_to_water) ** (1 / 3)  # number density according to water

        print('Build coordinates using Packmol: %s molecules ...' % self.n_mol_list)
        if self.procedure == Procedure.NPT_BINARY_SLAB:
            self.packmol.build_box(pdb_list, self.n_mol_list, 'init.pdb', length=length - 2, slab=True, silent=True)
        else:
            self.packmol.build_box(pdb_list, self.n_mol_list, 'init.pdb', length=length - 2, silent=True)

        print('Create box using DFF ...')
        if self.procedure == Procedure.NVT_SLAB:
            self.dff.build_bulk_after_packmol(pdb_list, self.n_mol_list, self.msd, pdb_corr='init.pdb',
                                              size=[length, length, length * 4])
        else:
            self.dff.build_bulk_after_packmol(pdb_list, self.n_mol_list, self.msd, pdb_corr='init.pdb',
                                              length=length)
