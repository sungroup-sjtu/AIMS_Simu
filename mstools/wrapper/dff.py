import os
import subprocess
import sys
from subprocess import PIPE, Popen

from mstools.errors import DffError


class DFF:
    TEMPLATE_DIR = os.path.abspath(os.path.dirname(__file__) + os.sep + '../template/dff/')
    '''
    wrappers for DFF
    '''
    pass

    def __init__(self, dff_root):
        self.DFF_ROOT = os.path.abspath(dff_root)
        if sys.platform == 'darwin':
            self.DFF_BIN_DIR = os.path.join(self.DFF_ROOT, 'bin32m')
        elif sys.platform.startswith('linux'):
            self.DFF_BIN_DIR = os.path.join(self.DFF_ROOT, 'bin32x')
        elif sys.platform.startswith('win'):
            self.DFF_BIN_DIR = os.path.join(self.DFF_ROOT, 'bin32w')
        else:
            raise DffError('Unsupported platform')

        self.DFFJOB_BIN = os.path.join(self.DFF_BIN_DIR, 'dffjob.exe')
        self.DFFEXP_BIN = os.path.join(self.DFF_BIN_DIR, 'dffexp.exe')

    def convert_model_to_msd(self, model, msd_out):
        pass

    def checkout(self, models: [str], db=None, table='TEAM_LS', ppf_out=None):
        if db == None:
            db = os.path.join(self.DFF_ROOT, 'database/TEAMFF.dffdb')
        if ppf_out == None:
            ppf_out = table + '.ppf'
        model_path = []
        for model in models:
            model_path.append(os.path.abspath(model))
        db = os.path.abspath(db)
        ppf_out = os.path.abspath(ppf_out)
        dfi = open(os.path.join(DFF.TEMPLATE_DIR, 't_checkout.dfi')).read()
        dfi = dfi.replace('%DATABASE%', db).replace('%TABLE%', table).replace('%MODELS%', '\n'.join(model_path)) \
            .replace('%OUTPUT%', ppf_out).replace('%LOG%', os.path.abspath('checkout.dfo'))
        with open('checkout.dfi', 'w') as f:
            f.write(dfi)
        sp = subprocess.Popen([self.DFFJOB_BIN, 'checkout'], stdin=PIPE, stdout=PIPE, stderr=PIPE)
        out, err = sp.communicate()
        if err.decode() != '':
            raise DffError('Checkout failed: %s' % err.decode())

    def set_charge(self, models: [str], ppf):
        model_path = []
        for model in models:
            model_path.append(os.path.abspath(model))
        ppf = os.path.abspath(ppf)
        log = os.path.abspath('setcharge.dfo')
        dfi = open(os.path.join(DFF.TEMPLATE_DIR, 't_set_charge.dfi')).read()
        dfi = dfi.replace('%MODELS%', '\n'.join(model_path)).replace('%PPF%', ppf).replace('%LOG%', log)
        with open('setcharge.dfi', 'w') as f:
            f.write(dfi)
        sp = subprocess.Popen([self.DFFJOB_BIN, 'setcharge'], stdin=PIPE, stdout=PIPE, stderr=PIPE)
        out, err = sp.communicate()
        if err.decode() != '':
            raise DffError('Set charge failed: %s' % err.decode())

    def export_lammps(self, model, ppf, data_out='data', lmp_out='in.lmp'):
        model = os.path.abspath(model)
        ppf = os.path.abspath(ppf)
        data_out = os.path.abspath(data_out)
        lmp_out = os.path.abspath(lmp_out)

        fftype = self.get_ff_type_from_ppf(ppf)

        dfi = open(os.path.join(DFF.TEMPLATE_DIR, 't_export_lammps.dfi')).read()
        dfi = dfi.replace('%ROOT%', self.DFF_ROOT).replace('%MODEL%', model).replace('%PPF%', ppf) \
            .replace('%TEMPLATE%', 'LAMMPS.Opt').replace('%FFTYPE%', fftype) \
            .replace('%DATAFILE%', data_out).replace('%INFILE%', lmp_out)
        with open('export_lammps.dfi', 'w') as f:
            f.write(dfi)
        sp = subprocess.Popen([self.DFFEXP_BIN, 'export_lammps'], stdin=PIPE, stdout=PIPE, stderr=PIPE)
        out, err = sp.communicate()
        if err.decode() != '':
            raise DffError('Export failed: %s' % err.decode())

    def export_gmx(self, model, ppf, gro_out='conf.gro', top_out='topol.top', mdp_out='grompp.mdp'):
        model = os.path.abspath(model)
        ppf = os.path.abspath(ppf)
        gro_out = os.path.abspath(gro_out)
        top_out = os.path.abspath(top_out)
        mdp_out = os.path.abspath(mdp_out)
        itp_out = top_out[:-4] + '.itp'

        fftype = self.get_ff_type_from_ppf(ppf)

        dfi = open(os.path.join(DFF.TEMPLATE_DIR, 't_export_gmx.dfi')).read()
        dfi = dfi.replace('%ROOT%', self.DFF_ROOT).replace('%MODEL%', model).replace('%PPF%', ppf) \
            .replace('%TEMPLATE%', 'GROMACS.Opt').replace('%FFTYPE%', fftype) \
            .replace('%GROFILE%', gro_out).replace('%TOPFILE%', top_out).replace('%MDPFILE%', mdp_out) \
            .replace('%ITPFILE%', itp_out)
        with open('export_gmx.dfi', 'w') as f:
            f.write(dfi)
        sp = Popen([self.DFFEXP_BIN, 'export_gmx'], stdin=PIPE, stdout=PIPE, stderr=PIPE)
        out, err = sp.communicate()
        if err.decode() != '':
            raise DffError('Export failed: %s' % err.decode())

    def build_box_after_packmol(self, models: [str], numbers: [int], msd_out, mol_corr,
                                size: [float] = None, length: float = None):
        if size != None:
            if len(size) != 3:
                raise DffError('Invalid box size')
            else:
                box = size
        elif length != None:
            box = [length, length, length]
        else:
            raise DffError('Box size needed')

        model_list_str = ''
        for model in models:
            model = os.path.abspath(model)
            model_list_str += '  MOL=%s' % model
        number_list_str = ' '.join(map(str, numbers))

        msd_out = os.path.abspath(msd_out)
        mol_corr = os.path.abspath(mol_corr)
        dfi = open(os.path.join(DFF.TEMPLATE_DIR, 't_packmol_multiple.dfi')).read()
        dfi = dfi.replace('%MODEL_LIST%', model_list_str).replace('%NUMBER_LIST%', number_list_str) \
            .replace('%OUT%', msd_out).replace('%CORRMOL%', mol_corr).replace('%PBC%', ' '.join(map(str, box)))
        with open('build.dfi', 'w') as f:
            f.write(dfi)
        sp = Popen([self.DFFJOB_BIN, 'build'], stdin=PIPE, stdout=PIPE, stderr=PIPE)
        out, err = sp.communicate()
        if err.decode() != '':
            raise DffError('Build failed: %s' % err.decode())

    @staticmethod
    def get_ff_type_from_ppf(ppf):
        with open(ppf) as f_ppf:
            for line in f_ppf:
                if line.startswith('#PROTOCOL'):
                    fftype = line.split('=')[-1].strip()
                    return fftype
        raise DffError('Unknown FF Type: %s' % ppf)
