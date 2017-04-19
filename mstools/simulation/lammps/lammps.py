import math
import os
import shutil

from ..simulation import Simulation
from ...errors import LammpsError
from ...jobmanager import Local
from ...utils import create_pdb_from_smiles
from ...wrapper import Packmol, DFF, Lammps


class LammpsSimulation(Simulation):
    def __init__(self, packmol_bin=None, dff_root=None, lmp_bin=None, jobmanager=Local()):
        super().__init__(packmol_bin=packmol_bin, dff_root=dff_root, jobmanager=jobmanager)
        self.LMP_BIN = lmp_bin

    def build(self, smiles, n_atoms=3000, density=1.0, ff='TEAM_LS',
              data_out='data', in_lmp='em.lmp',
              minimize=False):
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

    def prepare(self):
        shutil.copy('build/em.data', 'em.data')
        special, bond, angle, dihedral, improper = Lammps.get_intra_style_from_lmp('build/em.lmp')
        Lammps.prepare_lmp_from_template('t_npt.lmp', 'in.lmp', 'em.data', self.t, self.p / Unit.bar, int(1E3),
                                         self.n_mol, special, bond, angle, dihedral, improper)
