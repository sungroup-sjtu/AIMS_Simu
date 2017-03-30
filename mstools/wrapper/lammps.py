import os, subprocess, random
from subprocess import PIPE
from typing import List


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
            subprocess.Popen([self.LMP_BIN, '-i', in_lmp], stdout=PIPE).communicate()
        else:
            subprocess.Popen([self.LMP_BIN, '-i', in_lmp]).communicate()

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

    def rst_to_data(self, rst, data_out):
        subprocess.Popen([self.LMP_BIN, '-r', rst, data_out], stdout=PIPE).communicate()

    def get_box(self, data) -> List[float]:
        box = [0, 0, 0]
        with open(data) as f_data:
            for line in f_data:
                words = line.strip().split()
                if line.strip().endswith('xhi'):
                    box[0] = float(words[1]) - float(words[0])
                if line.strip().endswith('yhi'):
                    box[1] = float(words[1]) - float(words[0])
                if line.strip().endswith('zhi'):
                    box[2] = float(words[1]) - float(words[0])
                    break
        return box

    def scale_box(self, data, data_out, scale: List(float), remap=False):
        with open(data) as f_data:
            with open(data_out) as f_out:
                ATOMS_START = False
                for line in f_data:
                    if line.strip() == '':
                        f_out.write(line)
                        continue
                    if line.startswith('Atoms'):
                        ATOMS_START = True
                        f_out.write(line)
                        continue
                    if ATOMS_START and (not line.strip()[0].isdigit()):
                        ATOMS_START = False
                        f_out.write(line)
                        continue



