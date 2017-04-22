import os
import shutil
import subprocess
from subprocess import Popen, PIPE

from mstools.errors import GmxError


class GMX:
    TEMPLATE_DIR = os.path.abspath(os.path.dirname(__file__) + os.sep + '../template/gmx/')
    '''
    wrappers for GROMACS
    '''
    pass

    def __init__(self, gmx_bin):
        self.GMX_BIN = gmx_bin

    def grompp(self, mdp='grompp.mdp', gro='conf.gro', top='topol.top', tpr_out='md.tpr',
               maxwarn=3, silent=False, get_cmd=False):
        cmd = '%s grompp -f %s -c %s -p %s -o %s -maxwarn %i' % (self.GMX_BIN, mdp, gro, top, tpr_out, maxwarn)
        if get_cmd:
            return cmd
        else:
            (stdout, stderr) = (PIPE, PIPE) if silent else (None, None)
            sp = Popen(cmd.split(), stdout=stdout, stderr=stderr)
            sp.communicate()

    def mdrun(self, name='md', nprocs=1, rerun: str = None, silent=False, get_cmd=False):
        cmd = '%s mdrun -deffnm %s' % (self.GMX_BIN, name)
        if nprocs > 1:
            cmd = 'mpirun -np %i %s' % (nprocs, cmd)

        if rerun != None:
            cmd = '%s -rerun %s' % (cmd, rerun)

        if get_cmd:
            return cmd
        else:
            (stdout, stderr) = (PIPE, PIPE) if silent else (None, None)
            sp = Popen(cmd.split(), stdout=stdout, stderr=stderr)
            sp.communicate()

    def minimize(self, gro, top, nprocs=1, silent=False, name='em', pbc=True):
        if pbc:
            self.prepare_mdp_from_template('t_em.mdp')
        else:
            self.prepare_mdp_from_template('t_em_vacuum.mdp')

        self.grompp(gro=gro, top=top, tpr_out=name + '.tpr', silent=silent)
        self.mdrun(name=name, nprocs=nprocs, silent=silent)

    @staticmethod
    def prepare_mdp_from_template(template: str, mdp_out='grompp.mdp', T=298, P=1, nsteps=1000, dt=0.001,
                                  nstenergy=1000, nstvout=0,
                                  nstxtcout=10000, xtcgrps='System',
                                  genvel=True):

        if T == None or P == None or nsteps == None:
            raise GmxError('T, P, nsteps are required')

        genvel = 'yes' if genvel else 'no'

        template = os.path.join(GMX.TEMPLATE_DIR, template)
        if not os.path.exists(template):
            raise GmxError('mdp template not found')

        with open(template) as f_t:
            with open(mdp_out, 'w') as f_mdp:
                f_mdp.write(
                    f_t.read().replace('%T%', str(T)).replace('%P%', str(P)).replace('%nsteps%', str(int(nsteps))) \
                        .replace('%dt%', str(dt)) \
                        .replace('%nstenergy%', str(nstenergy)).replace('%nstvout%', str(nstvout)) \
                        .replace('%nstxtcout%', str(nstxtcout)).replace('%xtcgrps%', str(xtcgrps)) \
                        .replace('%genvel%', str(genvel)))

    def energy(self, edr, properties: [str], begin=0, get_cmd=False):
        property_str = '\\n'.join(properties)
        cmd = '%s energy -f %s -b %s' % (self.GMX_BIN, edr, str(begin))
        if get_cmd:
            cmd = 'echo "%s" | %s' % (property_str, cmd)
            return cmd
        else:
            sp = Popen(cmd.split(), stdout=PIPE, stdin=PIPE, stderr=PIPE)
            sp_out = sp.communicate(input=property_str.encode())[0]
            return sp_out

    def get_property(self, edr, property: str, begin=0) -> float:
        sp_out = self.energy(edr, properties=[property], begin=begin)

        for line in sp_out.decode().splitlines():
            if line.lower().startswith(property.lower()):
                return float(line.split()[1])
        raise GmxError('Invalid property')

    def get_box(self, edr, begin=0) -> [float]:
        sp = subprocess.Popen([self.GMX_BIN, 'energy', '-f', edr, '-b', str(begin)], stdout=PIPE, stdin=PIPE,
                              stderr=PIPE)
        sp_out = sp.communicate(input='Box-'.encode())[0]

        box = [0, 0, 0]
        for line in sp_out.decode().splitlines():
            if line.startswith('Box-X'):
                box[0] = float(line.split()[1])
            if line.startswith('Box-Y'):
                box[1] = float(line.split()[1])
            if line.startswith('Box-Z'):
                box[2] = float(line.split()[1])
        return box

    @staticmethod
    def scale_box(gro, new_gro, new_box: [float]):
        n = 0
        NAtoms = 0
        nresidue = []
        residue = []
        element = []
        natom = []
        xyz = []
        vxyz = []
        box = []
        if not os.path.exists(gro):
            gro = gro + '.gro'
        if not os.path.exists(gro):
            raise GmxError('gro not found')
        with open(gro) as f:
            for line in f:
                n += 1
                if n == 1:
                    continue
                if n == 2:
                    NAtoms = int(line.strip())
                    continue
                if n > 2 and n <= 2 + NAtoms:
                    nresidue.append(int(line[:5]))
                    residue.append(line[5:10])
                    element.append(line[10:15])
                    natom.append(int(line[15:20]))
                    x = float(line[20:28])
                    y = float(line[28:36])
                    z = float(line[36:44])
                    vx = float(line[44:52])
                    vy = float(line[52:60])
                    vz = float(line[60:68])
                    xyz.append([x, y, z])
                    vxyz.append([vx, vy, vz])
                if n == 3 + NAtoms:
                    box = [float(word) for word in line.strip().split()[:3]]
                    break

        scale = [new_box[i] / box[i] for i in range(3)]
        xyz = [[i[0] * scale[0], i[1] * scale[1], i[2] * scale[2]] for i in xyz]
        vxyz = [[i[0] * scale[0], i[1] * scale[1], i[2] * scale[2]] for i in vxyz]

        with open(new_gro, 'w') as f:
            f.write('Scaled box\n%i\n' % NAtoms)
            for i in range(NAtoms):
                f.write('%5i%5s%5s%5i%8.3f%8.3f%8.3f%8.4f%8.4f%8.4f\n'
                        % (nresidue[i], residue[i], element[i], natom[i],
                           xyz[i][0], xyz[i][1], xyz[i][2], vxyz[i][0], vxyz[i][1], vxyz[i][2]))
            f.write('%f %f %f\n' % (new_box[0], new_box[1], new_box[2]))

    @staticmethod
    def generate_top(itp, molecules, numbers):
        shutil.copy(itp, 'topol.top')
        with open('topol.top', 'a') as f:
            f.write('\n[system]\n[molecules]\n')
            for i, molecule in enumerate(molecules):
                f.write('%s %i\n' % (molecule, numbers[i]))

    def pdb2gro(self, pdb, gro_out, box: [float]):
        if len(box) != 3:
            raise GmxError('Invalid box')
        sp = Popen([self.GMX_BIN, 'editconf', '-f', pdb, '-o', gro_out, '-box', str(box[0]), str(box[1]), str(box[2])])
        sp.communicate()

    def velacc(self, trr, tpr=None, group=None, begin=0, xvg_out='velacc', silent=False):
        if tpr == None:
            tpr = trr
        if group == None:
            raise GmxError('No group specifed')

        (stdout, stderr) = (PIPE, PIPE) if silent else (None, None)
        sp = Popen([self.GMX_BIN, 'velacc', '-f', trr, '-s', tpr, '-o', xvg_out, '-b', str(begin), '-mol',
                    '-nonormalize'],
                   stdin=PIPE, stdout=stdout, stderr=stderr)
        sp.communicate(input=str(group).encode())

    def msd_com(self, xtc, tpr, resname, beginfit=-1, endfit=-1, xvg_out=None, silent=False):
        ndx = 'com-' + resname + '.ndx'
        GMX.select_com(tpr, resname, ndx_out=ndx)
        if xvg_out == None:
            xvg_out = 'msd-com-%s.xvg' % resname

        (stdout, stderr) = (PIPE, PIPE) if silent else (None, None)
        sp = Popen([self.GMX_BIN, 'msd', '-f', xtc, '-s', tpr, '-n', ndx, '-o', xvg_out, '-nomw',
                    '-beginfit', str(beginfit), '-endfit', str(endfit)], stdout=stdout, stderr=stderr)
        sp_out = sp.communicate()[0]

        for line in sp_out.decode().splitlines():
            if line.startswith('D['):
                return line
        raise GmxError('Error running gmx msd')

    def traj_com(self, xtc, tpr, xtc_out=None, silent=False):
        ndx = 'com.ndx'
        self.select_com(tpr, 'all', ndx_out=ndx)
        if xtc_out == None:
            xtc_out = 'com.xtc'

        (stdout, stderr) = (PIPE, PIPE) if silent else (None, None)
        sp = Popen([self.GMX_BIN, 'traj', '-s', tpr, '-f', xtc, '-oxt', xtc_out, '-mol', '-n', ndx],
                   stdout=stdout, stderr=stderr)
        sp.communicate()
        return xtc_out, ndx

    def select_com(self, tpr, resname, ndx_out='selection.ndx'):
        sp = Popen([self.GMX_BIN, 'select', '-s', tpr, '-on', ndx_out], stdin=PIPE, stdout=PIPE, stderr=PIPE)
        if resname == 'all':
            select_com_str = 'res_com of all'
        else:
            select_com_str = 'res_com of resname %s' % resname
        sp.communicate(input=select_com_str.encode())

    @staticmethod
    def generate_top_for_hvap(top, top_out):
        with open(top) as f:
            lines = f.read().splitlines()
        lines = [l for l in lines if not (l.startswith(';') or l == '')]

        f_out = open(top_out, 'w')

        PAIR_SECTION = False
        ATOMS_SECTION = False
        n_atoms = 0
        for line in lines:
            if line.find('[') != -1:
                ATOMS_SECTION = False
                PAIR_SECTION = False
            if line.find('[ atoms ]') != -1:
                ATOMS_SECTION = True
                n_atoms = 0
                f_out.write(line + '\n')
                continue
            if line.find('[ pairs ]') != -1:
                PAIR_SECTION = True
                f_out.write('[ exclusions ]\n')
                for i in range(1, n_atoms + 1):
                    other_atoms = list(range(1, n_atoms + 1))
                    other_atoms.remove(i)
                    exclusions = [i] + other_atoms
                    f_out.write(' '.join(list(map(str, exclusions))) + '\n')
            if ATOMS_SECTION:
                n_atoms += 1
            if not PAIR_SECTION:
                f_out.write(line + '\n')

        f_out.close()
