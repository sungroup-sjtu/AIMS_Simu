import pybel
from sqlalchemy import Column, ForeignKey, Integer, Text, String, Boolean, DateTime, and_, Float

from . import db

from functools import partial

NotNullColumn = partial(Column, nullable=False)


def is_exp_from_code(code):
    return str(code).find('1') != -1


class Group(db.Model):
    __bind_key__ = 'yaws'
    __tablename__ = 'group'
    id = NotNullColumn(Integer, primary_key=True)
    name = NotNullColumn(Text)
    smarts = NotNullColumn(Text)

    molecules = db.relationship('YawsMolecule', secondary='molecule_group')

    def __repr__(self):
        return '<Group: %s %s>' % (self.name, self.smarts)


class YawsMolecule(db.Model):
    __bind_key__ = 'yaws'
    __tablename__ = 'yaws_molecule'
    id = NotNullColumn(Integer, primary_key=True)
    formula = NotNullColumn(Text)
    name = NotNullColumn(Text)
    cas = Column(Text, nullable=True)
    weight = Column(Float, nullable=True)
    Tfus = Column(Float, nullable=True)
    Tfus_code = Column(Text, nullable=True)
    Tvap = Column(Float, nullable=True)
    Tvap_code = Column(Text, nullable=True)
    cid = Column(Integer, nullable=True)
    iupac = Column(Text, nullable=True)
    smiles = Column(Text, nullable=True)
    isomeric_smiles = Column(Text, nullable=True)
    category = Column(Text, nullable=True)

    groups = db.relationship('Group', secondary='molecule_group')

    def __repr__(self):
        return '<Molecule: %s %s %s>' % (self.formula, self.isomeric_smiles, self.name)

    @property
    def groups_str(self):
        groups = [group.name for group in self.groups]
        return ', '.join(groups)

    @property
    def Tc(self):
        critical = YawsCritical.query.filter_by(cas=self.cas).first()
        if critical is None:
            return None
        return critical.Tc

    @property
    def n_heavy_atom(self):
        return pybel.readstring('smi', self.smiles).OBMol.NumHvyAtoms()

    def get_Tfus(self):
        return self.Tfus

    @property
    def Tfus_is_exp(self):
        return is_exp_from_code(self.Tfus_code)

    def get_Tvap(self):
        return self.Tvap

    @property
    def Tvap_is_exp(self):
        return is_exp_from_code(self.Tvap_code)

    def get_Tc(self):
        return self.Tc

    @property
    def Tc_is_exp(self):
        critical = YawsCritical.query.filter_by(cas=self.cas).first()
        if critical is not None:
            return is_exp_from_code(critical.Tc_code)
        return False

    def get_Pvap(self, T=298):
        row = YawsPvap.query.filter_by(cas=self.cas).first()
        if row is not None:
            return row.get_value_at_T(T)
        return None

    @property
    def Pvap_is_exp(self):
        row = YawsPvap.query.filter_by(cas=self.cas).first()
        if row is not None:
            return is_exp_from_code(row.code)
        return False

    def get_compressibility_T(self) -> (int, float):
        row = YawsCompressibility.query.filter_by(cas=self.cas).first()
        if row is not None:
            T = row.get_T()
            val = row.get_value_at_T(T)
            return T, val
        return None, None

    @property
    def compressibility_is_exp(self):
        row = YawsCompressibility.query.filter_by(cas=self.cas).first()
        if row is not None:
            return is_exp_from_code(row.code)
        return False

    def get_density(self, T=298):
        row = YawsDensity.query.filter_by(cas=self.cas).first()
        if row is not None:
            return row.get_value_at_T(T)
        return None

    @property
    def density_is_exp(self):
        row = YawsDensity.query.filter_by(cas=self.cas).first()
        if row is not None:
            return is_exp_from_code(row.code)
        return False

    def get_Hvap(self, T=298):
        row = YawsHvap.query.filter_by(cas=self.cas).first()
        if row is not None:
            return row.get_value_at_T(T)
        return None

    @property
    def Hvap_is_exp(self):
        row = YawsHvap.query.filter_by(cas=self.cas).first()
        if row is not None:
            return is_exp_from_code(row.code)
        return False

    def get_expansion(self, T=298):
        row = YawsExpansion.query.filter_by(cas=self.cas).first()
        if row is not None:
            return row.get_value_at_T(T)
        return None

    @property
    def expansion_is_exp(self):
        row = YawsExpansion.query.filter_by(cas=self.cas).first()
        if row is not None:
            return is_exp_from_code(row.code)
        return False

    def get_Cp(self, T=298):
        row = YawsCp.query.filter_by(cas=self.cas).first()
        if row is not None:
            return row.get_value_at_T(T)
        return None

    @property
    def Cp_is_exp(self):
        row = YawsCp.query.filter_by(cas=self.cas).first()
        if row is not None:
            return is_exp_from_code(row.code)
        return False

    def get_surface_tension(self, T=298):
        row = YawsST.query.filter_by(cas=self.cas).first()
        if row is not None:
            return row.get_value_at_T(T)
        return None

    @property
    def surface_tension_is_exp(self):
        row = YawsST.query.filter_by(cas=self.cas).first()
        if row is not None:
            return is_exp_from_code(row.code)
        return False

    def get_viscosity(self, T=298):
        row = YawsViscosity.query.filter_by(cas=self.cas).first()
        if row is not None:
            return row.get_value_at_T(T)
        return None

    @property
    def viscosity_is_exp(self):
        row = YawsViscosity.query.filter_by(cas=self.cas).first()
        if row is not None:
            return is_exp_from_code(row.code)
        return False

    def get_thermal_conductivity(self, T=298):
        row = YawsTC.query.filter_by(cas=self.cas).first()
        if row is not None:
            return row.get_value_at_T(T)
        return None

    @property
    def thermal_conductivity_is_exp(self):
        row = YawsTC.query.filter_by(cas=self.cas).first()
        if row is not None:
            return is_exp_from_code(row.code)
        return False


class MoleculeGroup(db.Model):
    __bind_key__ = 'yaws'
    __tablename__ = 'molecule_group'
    molecule_id = NotNullColumn(Integer, ForeignKey(YawsMolecule.id), primary_key=True)
    group_id = NotNullColumn(Integer, ForeignKey(Group.id), primary_key=True)


class YawsCompressibility(db.Model):
    __bind_key__ = 'yaws'
    __tablename__ = 'yaws_compressibility'
    id = NotNullColumn(Integer, primary_key=True)
    molecule_id = Column(Integer, ForeignKey(YawsMolecule.id), nullable=True)
    formula = NotNullColumn(Text)
    name = NotNullColumn(Text)
    cas = Column(Text)
    weight = NotNullColumn(Float)
    T = NotNullColumn(Float)
    value = NotNullColumn(Float)
    code = NotNullColumn(Text)

    def get_T(self):
        return self.T + 273.15

    def get_value_at_T(self, T):
        '''
        J/mol/K
        :param T:
        :return:
        '''
        if T < self.T + 273.15 - 1 or T > self.T + 273.15 + 1:
            return None
        return self.value


class YawsCp(db.Model):
    __bind_key__ = 'yaws'
    __tablename__ = 'yaws_Cp'
    id = NotNullColumn(Integer, primary_key=True)
    molecule_id = Column(Integer, ForeignKey(YawsMolecule.id), nullable=True)
    formula = NotNullColumn(Text)
    name = NotNullColumn(Text)
    cas = Column(Text)
    weight = NotNullColumn(Float)
    A = NotNullColumn(Float)
    B = NotNullColumn(Float)
    C = NotNullColumn(Float)
    D = NotNullColumn(Float)
    E = NotNullColumn(Float)
    Tmin = NotNullColumn(Float)
    Tmax = NotNullColumn(Float)
    code = NotNullColumn(Text)

    def get_value_at_T(self, T):
        '''
        J/mol/K
        :param T:
        :return:
        '''
        if T < self.Tmin - 1 or T > self.Tmax:
            return None
        value = self.A + self.B * T + self.C * T ** 2 + self.D * T ** 3 + self.E * T ** 4
        return value


class YawsCritical(db.Model):
    __bind_key__ = 'yaws'
    __tablename__ = 'yaws_critical'
    No = NotNullColumn(Integer, primary_key=True)
    formula = Column(Text)
    name = Column(Text)
    cas = Column(Text)
    weight = Column(Float)
    Tc = Column(Float)
    Tc_code = Column(Text)
    Pc = Column(Float)
    Pc_code = Column(Text)
    density = Column(Float)
    density_code = Column(Text)


class YawsDensity(db.Model):
    __bind_key__ = 'yaws'
    __tablename__ = 'yaws_density'
    id = NotNullColumn(Integer, primary_key=True)
    molecule_id = Column(Integer, ForeignKey(YawsMolecule.id), nullable=True)
    formula = NotNullColumn(Text)
    name = NotNullColumn(Text)
    cas = Column(Text)
    weight = NotNullColumn(Float)
    A = NotNullColumn(Float)
    B = NotNullColumn(Float)
    C = NotNullColumn(Float)
    n = NotNullColumn(Float)
    Tmin = NotNullColumn(Float)
    Tmax = NotNullColumn(Float)
    code = NotNullColumn(Text)

    def get_value_at_T(self, T):
        '''
        g/mL
        :param T:
        :return:
        '''
        if T < self.Tmin - 1 or T > self.Tmax:
            return None
        value = self.A * self.B ** (-(1 - T / self.C) ** self.n)
        return value


class YawsExpansion(db.Model):
    __bind_key__ = 'yaws'
    __tablename__ = 'yaws_expansion'
    id = NotNullColumn(Integer, primary_key=True)
    molecule_id = Column(Integer, ForeignKey(YawsMolecule.id), nullable=True)
    formula = NotNullColumn(Text)
    name = NotNullColumn(Text)
    cas = Column(Text)
    weight = NotNullColumn(Float)
    A = NotNullColumn(Float)
    B = NotNullColumn(Float)
    m = NotNullColumn(Float)
    Tmin = NotNullColumn(Float)
    Tmax = NotNullColumn(Float)
    code = NotNullColumn(Text)

    def get_value_at_T(self, T):
        '''
        /K
        :param T:
        :return:
        '''
        if T < self.Tmin - 1 or T > self.Tmax:
            return None
        value = self.A * (1 - T / self.B) ** self.m
        return value


class YawsHvap(db.Model):
    __bind_key__ = 'yaws'
    __tablename__ = 'yaws_Hvap'
    id = NotNullColumn(Integer, primary_key=True)
    molecule_id = Column(Integer, ForeignKey(YawsMolecule.id), nullable=True)
    formula = NotNullColumn(Text)
    name = NotNullColumn(Text)
    cas = Column(Text)
    weight = NotNullColumn(Float)
    A = NotNullColumn(Float)
    B = NotNullColumn(Float)
    n = NotNullColumn(Float)
    Tmin = NotNullColumn(Float)
    Tmax = NotNullColumn(Float)
    code = NotNullColumn(Text)

    def get_value_at_T(self, T):
        '''
        kJ/mol
        :param T:
        :return:
        '''
        if T < self.Tmin - 1 or T > self.Tmax:
            return None
        value = self.A * (1 - T / self.B) ** self.n
        return value


class YawsKow(db.Model):
    __bind_key__ = 'yaws'
    __tablename__ = 'yaws_Kow'
    id = NotNullColumn(Integer, primary_key=True)
    molecule_id = Column(Integer, ForeignKey(YawsMolecule.id), nullable=True)
    formula = NotNullColumn(Text)
    name = NotNullColumn(Text)
    cas = Column(Text)
    weight = NotNullColumn(Float)
    value = NotNullColumn(Float)
    code = NotNullColumn(Text)


class YawsPvap(db.Model):
    __bind_key__ = 'yaws'
    __tablename__ = 'yaws_Pvap'
    id = NotNullColumn(Integer, primary_key=True)
    molecule_id = Column(Integer, ForeignKey(YawsMolecule.id), nullable=True)
    formula = NotNullColumn(Text)
    name = NotNullColumn(Text)
    cas = Column(Text)
    weight = NotNullColumn(Float)
    A = NotNullColumn(Float)
    B = NotNullColumn(Float)
    C = NotNullColumn(Float)
    Tmin = NotNullColumn(Float)
    Tmax = NotNullColumn(Float)
    code = NotNullColumn(Text)

    def get_value_at_T(self, T):
        '''
        pascal
        :param T:
        :return:
        '''
        if T < self.Tmin or T > self.Tmax:
            return None
        value = 10 ** (self.A - self.B / (T - 273.15 + self.C)) * 0.00133322368
        return value


class YawsST(db.Model):
    __bind_key__ = 'yaws'
    __tablename__ = 'yaws_ST'
    id = NotNullColumn(Integer, primary_key=True)
    molecule_id = Column(Integer, ForeignKey(YawsMolecule.id), nullable=True)
    formula = NotNullColumn(Text)
    name = NotNullColumn(Text)
    cas = Column(Text)
    weight = NotNullColumn(Float)
    A = NotNullColumn(Float)
    B = NotNullColumn(Float)
    n = NotNullColumn(Float)
    Tmin = NotNullColumn(Float)
    Tmax = NotNullColumn(Float)
    code = NotNullColumn(Text)

    def get_value_at_T(self, T):
        '''
        mN/m
        :param T:
        :return:
        '''
        if T < self.Tmin - 1 or T > self.Tmax:
            return None
        value = self.A * (1 - T / self.B) ** self.n
        return value


class YawsTC(db.Model):
    __bind_key__ = 'yaws'
    __tablename__ = 'yaws_TC'
    id = NotNullColumn(Integer, primary_key=True)
    molecule_id = Column(Integer, ForeignKey(YawsMolecule.id), nullable=True)
    formula = NotNullColumn(Text)
    name = NotNullColumn(Text)
    cas = Column(Text)
    weight = NotNullColumn(Float)
    A = NotNullColumn(Float)
    B = NotNullColumn(Float)
    C = NotNullColumn(Float)
    Tmin = NotNullColumn(Float)
    Tmax = NotNullColumn(Float)
    code = NotNullColumn(Text)

    def get_value_at_T(self, T):
        '''
        W/m/K
        :param T:
        :return:
        '''
        if T < self.Tmin - 1 or T > self.Tmax:
            return None
        value = self.A + self.B * T + self.C * T ** 2
        return value


class YawsViscosity(db.Model):
    __bind_key__ = 'yaws'
    __tablename__ = 'yaws_viscosity'
    id = NotNullColumn(Integer, primary_key=True)
    molecule_id = Column(Integer, ForeignKey(YawsMolecule.id), nullable=True)
    formula = NotNullColumn(Text)
    name = NotNullColumn(Text)
    cas = Column(Text)
    weight = NotNullColumn(Float)
    A = NotNullColumn(Float)
    B = NotNullColumn(Float)
    C = NotNullColumn(Float)
    D = NotNullColumn(Float)
    Tmin = NotNullColumn(Float)
    Tmax = NotNullColumn(Float)
    code = NotNullColumn(Text)

    def get_value_at_T(self, T):
        '''
        cP
        :param T:
        :return:
        '''
        if T < self.Tmin - 1 or T > self.Tmax:
            return None
        value = 10 ** (self.A + self.B / T + self.C * T + self.D * T ** 2)
        return value


class YawsVLE(db.Model):
    __bind_key__ = 'yaws'
    __tablename__ = 'yaws_VLE'
    No = NotNullColumn(Integer, primary_key=True)
    formula = Column(Text)
    name = Column(Text)
    cas = Column(Text)
    weight = Column(Float)
    phase = Column(Text)
    T = Column(Float)
    P = Column(Float)
    density = Column(Float)
    code = Column(Text)
