import os
import shutil

from mstools.errors import GmxError
from mstools.simulation import Simulation
from mstools.wrapper import GMX


class GmxSimulation(Simulation):
    def __init__(self, packmol_bin=None, dff_root=None, gmx_bin=None, procedure=None):
        super().__init__(packmol_bin=packmol_bin, dff_root=dff_root, procedure=procedure)
        self.gmx = GMX(gmx_bin=gmx_bin)

    def build(self, smiles, n_atoms=3000, ff='TEAM_LS',
              gro_out='conf.gro', top_out='topol.top', mdp_out='grompp.mdp',
              minimize=False):
        self.build_dff_box_from_smiles(smiles, n_atoms)
        print('Checkout force field: %s ...' % ff)
        self.dff.checkout(self.msd, table=ff)
        print('Export GROMACS files ...')
        self.dff.export_gmx(self.msd, ff + '.ppf', gro_out, top_out, mdp_out)

        if minimize:
            print('Energy minimizing ...')
            self.gmx.minimize(gro_out, top_out, name='em', silent=True)

            if os.path.exists('em.gro'):
                shutil.move('em.gro', gro_out)
            else:
                raise GmxError('Energy minimization failed')
