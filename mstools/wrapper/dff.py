import os, subprocess, sys
from subprocess import PIPE, STDOUT


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
            raise Exception('DFF: Unsupported platform')

        self.DFFJOB_BIN = os.path.join(self.DFF_BIN_DIR, 'dffjob.exe')
        self.DFFEXP_BIN = os.path.join(self.DFF_BIN_DIR, 'dffexp.exe')

    def convert_model_to_msd(self, model, msd_out):
        pass

    def checkout(self, model, db, table, ppf_out):
        model = os.path.abspath(model)
        db = os.path.abspath(db)
        ppf_out = os.path.abspath(ppf_out)
        dfi = open(os.path.join(DFF.TEMPLATE_DIR, 't_checkout.dfi')).read()
        dfi = dfi.replace('%DATABASE%', db).replace('%TABLE%', table).replace('%MODEL%', model) \
            .replace('%OUTPUT%', ppf_out).replace('%LOG%', os.path.abspath('checkout.dfo'))
        with open('checkout.dfi', 'w') as f:
            f.write(dfi)
        subprocess.Popen([self.DFFJOB_BIN, 'checkout'], stdin=PIPE, stdout=PIPE).communicate()

    def export_lammps(self, model, ppf, data_out, lmp_out):
        model = os.path.abspath(model)
        ppf = os.path.abspath(ppf)
        data_out = os.path.abspath(data_out)
        lmp_out = os.path.abspath(lmp_out)

        fftype = 'TEAM'
        with open(ppf) as f_ppf:
            for line in f_ppf:
                if line.startswith('#PROTOCOL'):
                    fftype = line.split('=')[-1].strip()
                    break

        dfi = open(os.path.join(DFF.TEMPLATE_DIR, 't_export_lammps.dfi')).read()
        dfi = dfi.replace('%ROOT%', self.DFF_ROOT).replace('%MODEL%', model).replace('%PPF%', ppf) \
            .replace('%TEMPLATE%', 'LAMMPS.Opt').replace('%FFTYPE%', fftype) \
            .replace('%DATAFILE%', data_out).replace('%INFILE%', lmp_out)
        with open('export_lammps.dfi', 'w') as f:
            f.write(dfi)
        subprocess.Popen([self.DFFEXP_BIN, 'export_lammps'], stdin=PIPE, stdout=PIPE).communicate()

    def export_gmx(self, model, ppf, gro_out, top_out, itp_out, mdp_out):
        pass

    def build_bulk_after_packmol(self, model, number, msd_out, pdb_corr, length):
        model = os.path.abspath(model)
        msd_out = os.path.abspath(msd_out)
        pdb_corr = os.path.abspath(pdb_corr)
        dfi = open(os.path.join(DFF.TEMPLATE_DIR, 't_packmol.dfi')).read()
        dfi = dfi.replace('%MODEL%', model).replace('%NUMBER%', str(number)).replace('%OUT%', msd_out) \
            .replace('%PDB%', pdb_corr).replace('%LENGTH%', str(length))
        with open('build.dfi', 'w') as f:
            f.write(dfi)
        subprocess.Popen([self.DFFJOB_BIN, 'build'], stdin=PIPE, stdout=PIPE).communicate()

    def checkout_TEAM_MS(self, model, ppf_out):
        db = os.path.join(self.DFF_ROOT, 'database/TEAMFF.dffdb')
        table = 'TEAM_MS'
        self.checkout(model, db, table, ppf_out)

    def checkout_TEAM_LS(self, model, ppf_out):
        db = os.path.join(self.DFF_ROOT, 'database/TEAMFF.dffdb')
        table = 'TEAM_LS'
        self.checkout(model, db, table, ppf_out)
