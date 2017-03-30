import os, subprocess, random
from subprocess import PIPE


class Lammps:
    TEMPLATE_DIR = os.path.abspath(os.path.dirname(__file__) + os.sep + '../template/lammps/')
    '''
    wrappers for Lammps
    '''
    pass

    def __init__(self, lmp_bin):
        self.LMP_BIN = lmp_bin

    def run(self, in_lmp, silent=False):
        if silent:
            subprocess.check_call([self.LMP_BIN, '-i', in_lmp], stdout=PIPE)
        else:
            subprocess.check_call([self.LMP_BIN, '-i', in_lmp])

    @staticmethod
    def prepare_lmp_from_template(template: str, lmp, datafile, T, P, nsteps, n_mol=1, specialbonds=None,
                                  bondstyle=None, anglestyle=None, dihedralstyle=None, improperstyle=None):

        template = os.path.join(Lammps.TEMPLATE_DIR, template)
        if not os.path.exists(template):
            raise Exception('mdp template not found')

        specialbonds = specialbonds or 'lj/coul 0 0 1'
        bondstyle = bondstyle or 'none'
        anglestyle = anglestyle or 'none'
        dihedralstyle = dihedralstyle or 'none'
        improperstyle = improperstyle or 'none'

        with open(template) as f_t:
            with open(lmp, 'w') as f_lmp:
                f_lmp.write(
                    f_t.read().replace('%T%', str(T)).replace('%P%', str(P)).replace('%STEPS%', str(int(nsteps)))
                    .replace('%NMOL%', str(n_mol))
                    .replace('%DATAFILE%', datafile).replace('%SPECIALBONDS%', specialbonds)
                    .replace('%BONDSTYLE%', bondstyle).replace('%ANGLESTYLE%', anglestyle)
                    .replace('%DIHEDRALSTYLE%', dihedralstyle).replace('%IMPROPERSTYLE%', improperstyle)
                    .replace('%RANDINT%', str(random.randint(1E7, 1E8))))

    @staticmethod
    def get_intra_style_from_lmp(lmp):
        special = None
        bond = None
        angle = None
        dihedral = None
        improper = None
        with open(lmp) as f_lmp:
            for line in f_lmp:
                if line.startswith('special_bonds'):
                    special = ' '.join(line.strip().split()[1:])
                elif line.startswith('bond_style'):
                    bond = ' '.join(line.strip().split()[1:])
                elif line.startswith('angle_style'):
                    angle = ' '.join(line.strip().split()[1:])
                elif line.startswith('dihedral_style'):
                    dihedral = ' '.join(line.strip().split()[1:])
                elif line.startswith('improper_style'):
                    improper = ' '.join(line.strip().split()[1:])
        return special, bond, angle, dihedral, improper
