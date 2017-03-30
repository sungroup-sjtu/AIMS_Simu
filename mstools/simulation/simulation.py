import os
import math
import shutil
from mstools.wrapper import Packmol, DFF, Lammps


class Simulation():
    def __init__(self, packmol_bin=None, dff_root=None, lmp_bin=None):
        self.PACKMOL_BIN = packmol_bin
        self.LMP_BIN = lmp_bin
        self.DFF_ROOT = dff_root

    def build(self):
        pass

    def prepare(self):
        pass

    def run(self):
        pass

    def analyze(self):
        pass

    def build_lammps_box_from_smiles(self, smiles, n_atoms, data, in_lmp, density=0.8, ff='TEAM_MS', minimize=False):
        import pybel
        try:
            py_mol = pybel.readstring('smi', smiles)
            py_mol.addh()
            py_mol.make3D()
            py_mol.write('pdb', 'mol.pdb', overwrite=True)
        except:
            raise Exception('Cannot create PDB from SMILES')

        n_atom_per_mol = len(py_mol.atoms)
        number = math.ceil(n_atoms / n_atom_per_mol)  # A total of 3000 atoms
        length = (10 / 6.022 * py_mol.molwt * number / density) ** (1 / 3)

        self.n_mol = number
        packmol = Packmol(packmol_bin=self.PACKMOL_BIN)
        print('Build coordinates using packmol...')
        packmol.build_box(['mol.pdb'], [number], 'init.pdb', length=length - 2)
        dff = DFF(dff_root=self.DFF_ROOT)
        print('Create box using DFF...')
        dff.build_bulk_after_packmol('mol.pdb', number, 'init.msd', pdb_corr='init.pdb', length=length)
        print('Checkout force field...')
        if ff == 'TEAM_MS':
            dff.checkout_TEAM_MS('init.msd', 'TEAM_MS.ppf')
        elif ff == 'TEAM_LS':
            dff.checkout_TEAM_LS('init.msd', 'TEAM_LS.ppf')
        else:
            raise Exception('Unknown force field: %s' % ff)
        print('Export lammps data file...')
        dff.export_lammps('init.msd', 'TEAM_MS.ppf', data, in_lmp)

        if minimize:
            lammps = Lammps(self.LMP_BIN)
            print('Energy minimizing...')
            lammps.run(in_lmp, silent=True)

            if os.path.exists('em.data'):
                shutil.copy('em.data', data)
            else:
                raise Exception('Energy minimization failed')
