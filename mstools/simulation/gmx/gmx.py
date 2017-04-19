import os
import shutil

from ..procedure import Procedure
from ..simulation import Simulation
from ...errors import GmxError
from ...jobmanager import Local
from ...unit import Unit
from ...wrapper import GMX


class GmxSimulation(Simulation):
    def __init__(self, packmol_bin=None, dff_root=None, gmx_bin=None, jobmanager=Local()):
        super().__init__(packmol_bin=packmol_bin, dff_root=dff_root, jobmanager=jobmanager)
        self.gmx = GMX(gmx_bin=gmx_bin)

    def build(self, smiles_list: [str], n_atoms=3000, ff='TEAM_LS',
              gro_out='conf.gro', top_out='topol.top', mdp_out='grompp.mdp',
              minimize=False):
        self.build_dff_box_from_smiles(smiles_list, n_atoms)
        print('Checkout force field: %s ...' % ff)
        self.dff.checkout(self.msd, table=ff)
        print('Export GROMACS files ...')
        self.dff.export_gmx(self.msd, ff + '.ppf', gro_out, top_out, mdp_out)

        if minimize:
            print('Energy minimize ...')
            self.gmx.minimize(gro_out, top_out, name='em', silent=True)

            if os.path.exists('em.gro'):
                shutil.move('em.gro', gro_out)
            else:
                raise GmxError('Energy minimization failed')

    def prepare(self, model_dir='.', gro='conf.gro', top='topol.top', T=None, P=None, nproc=1, jobname=None):
        if os.path.abspath(model_dir) != os.getcwd():
            shutil.copy(os.path.join(model_dir, gro), gro)
            shutil.copy(os.path.join(model_dir, top), top)
            for f in os.listdir(model_dir):
                if f.endswith('.itp'):
                    shutil.copy(os.path.join(model_dir, f), '.')

        commands = []
        if self.procedure in (Procedure.NPT, Procedure.NPT_BINARY_SLAB,):
            self.gmx.prepare_mdp_from_template('t_npt.mdp', T=T, P=P/Unit.bar, nsteps=int(1E6))
        elif self.procedure in (Procedure.NVT_SLAB,):
            self.gmx.prepare_mdp_from_template('t_nvt.mdp', T=T, nsteps=int(1E6))
        cmd = self.gmx.grompp(gro=gro, top=top, tpr_out=self.procedure, get_cmd=True)
        commands.append(cmd)
        cmd = self.gmx.mdrun(name=self.procedure, nprocs=nproc, get_cmd=True)
        commands.append(cmd)
        self.jobmanager.generate_sh(os.getcwd(), commands, name=jobname or self.procedure)

    def analyze(self):
        import panedr
        df = panedr.edr_to_df(self.procedure + '.edr')
        pass
