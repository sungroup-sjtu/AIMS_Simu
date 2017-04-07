import shutil

from mstools.errors import *
from mstools.utils import *
from mstools.wrapper import *


class Simulation():
    def __init__(self, packmol_bin=None, dff_root=None, lmp_bin=None, gmx_bin=None):
        self.PACKMOL_BIN = packmol_bin
        self.DFF_ROOT = dff_root
        self.LMP_BIN = lmp_bin
        self.GMX_BIN = gmx_bin

    def build(self):
        pass

    def prepare(self):
        pass

    def run(self):
        pass

    def analyze(self):
        pass

    def build_lammps_box_from_smiles(self, smiles, n_atoms=3000, density=0.8,
                                     data_out='data', in_lmp='em.lmp', ff='TEAM_LS', minimize=False):
        py_mol = create_pdb_from_smiles(smiles, 'mol.pdb')
        n_atom_per_mol = len(py_mol.atoms)
        number = math.ceil(n_atoms / n_atom_per_mol)  # A total of 3000 atoms
        length = (10 / 6.022 * py_mol.molwt * number / density) ** (1 / 3)

        self.n_mol = number
        packmol = Packmol(packmol_bin=self.PACKMOL_BIN)
        print('Build coordinates using packmol: ~ %i atoms ...' % n_atoms)
        packmol.build_box(['mol.pdb'], [number], 'init.pdb', length=length - 2)
        dff = DFF(dff_root=self.DFF_ROOT)
        print('Create box using DFF...')
        dff.build_bulk_after_packmol('mol.pdb', number, 'init.msd', pdb_corr='init.pdb', length=length)
        print('Checkout force field: %s ...' % ff)
        dff.checkout('init.msd', table=ff)
        print('Export lammps files...')
        dff.export_lammps('init.msd', ff + '.ppf', data_out, in_lmp)

        if minimize:
            lammps = Lammps(self.LMP_BIN)
            print('Energy minimizing...')
            lammps.run(in_lmp, silent=True)

            if os.path.exists('em.data'):
                shutil.move('em.data', data_out)
            else:
                raise LammpsError('Energy minimization failed')

    def build_gmx_box_from_smiles(self, smiles, n_atoms=3000, density=0.8,
                                  gro_out='conf.gro', top_out='topol.top', mdp_out='grompp.mdp',
                                  ff='TEAM_LS', minimize=False):
        py_mol = create_pdb_from_smiles(smiles, 'mol.pdb')
        n_atom_per_mol = len(py_mol.atoms)
        number = math.ceil(n_atoms / n_atom_per_mol)  # A total of 3000 atoms
        length = (10 / 6.022 * py_mol.molwt * number / density) ** (1 / 3)

        self.n_mol = number
        packmol = Packmol(packmol_bin=self.PACKMOL_BIN)
        print('Build coordinates using packmol: ~ %i atoms ...' % n_atoms)
        packmol.build_box(['mol.pdb'], [number], 'init.pdb', length=length - 2)
        dff = DFF(dff_root=self.DFF_ROOT)
        print('Create box using DFF...')
        dff.build_bulk_after_packmol('mol.pdb', number, 'init.msd', pdb_corr='init.pdb', length=length)
        print('Checkout force field: %s ...' % ff)
        dff.checkout('init.msd', table=ff)
        print('Export GROMACS files...')
        dff.export_gmx('init.msd', ff + '.ppf', gro_out, top_out, mdp_out)

        if minimize:
            gmx = GMX(gmx_bin=self.GMX_BIN)
            print('Energy minimizing...')
            gmx.grompp(mdp_out, gro_out, top_out, tpr_out='em.tpr', silent=True)
            gmx.mdrun('em', silent=True)

            if os.path.exists('em.gro'):
                shutil.move('em.gro', gro_out)
            else:
                raise GmxError('Energy minimization failed')
