from . import db
from sqlalchemy import Column, Integer, Float, Text, String, ForeignKey
from sqlalchemy.orm import relationship

from functools import partial

NotNullColumn = partial(Column, nullable=False)


class NistMolecule(db.Model):
    __bind_key__ = 'nist'
    __tablename__ = 'nist_molecule'
    id = Column(Integer, primary_key=True)
    content_id = Column(String(100), unique=True)
    knovel_id = Column(Integer)
    cas = Column(String(100))  # CAS could be duplicated for isomers
    formula = Column(String(100))
    name = Column(Text)
    smiles = Column(Text)
    inchi = Column(Text)
    weight = Column(Float)
    tc = Column(Float)  # K
    pc = Column(Float)  # kPa
    dc = Column(Float)  # kg/m^3
    tb = Column(Float)  # K
    tt = Column(Float)  # K
    hfus = Column(Float)  # kJ/mol
    tc_u = Column(Float)  # K
    pc_u = Column(Float)  # kPa
    dc_u = Column(Float)  # kg/m^3
    tb_u = Column(Float)  # K
    tt_u = Column(Float)  # K
    hfus_u = Column(Float)  # kJ/mol
    remark = Column(Text)
    n_heavy = Column(Integer)

    datas = relationship('NistData', lazy='dynamic')
    splines = relationship('NistSpline', lazy='dynamic')

    def __repr__(self):
        return '<NistMolecule: %s %s %s>' % (self.formula, self.smiles, self.name)

    def update_n_heavy(self):
        import pybel
        m = pybel.readstring('smi', self.smiles)
        self.n_heavy = m.OBMol.NumHvyAtoms()


class NistProperty(db.Model):
    __bind_key__ = 'nist'
    __tablename__ = 'nist_property'
    id = Column(Integer, primary_key=True)
    name = Column(Text)
    property_id = Column(String(100))
    phase_id = Column(String(100))

    datas = relationship('NistData', lazy='dynamic')
    splines = relationship('NistSpline', lazy='dynamic')


class NistData(db.Model):
    __bind_key__ = 'nist'
    __tablename__ = 'nist_data'
    id = Column(Integer, primary_key=True)
    molecule_id = Column(Integer, ForeignKey(NistMolecule.id))
    property_id = Column(Integer, ForeignKey(NistProperty.id))
    t = Column(Float)  # K
    p = Column(Float)  # kPa
    value = Column(Float)
    uncertainty = Column(Float)

    molecule = relationship(NistMolecule)
    property = relationship(NistProperty)


class NistSpline(db.Model):
    __bind_key__ = 'nist'
    __tablename__ = 'nist_spline'
    id = Column(Integer, primary_key=True)
    molecule_id = Column(Integer, ForeignKey(NistMolecule.id))
    property_id = Column(Integer, ForeignKey(NistProperty.id))
    t_min = Column(Float)
    t_max = Column(Float)
    coef_v = Column(Text)  # value
    coef_u = Column(Text)  # uncertainty

    molecule = relationship(NistMolecule)
    property = relationship(NistProperty)

    def get_data(self, T):
        if T < self.t_min or T > self.t_max:
            return None, None

        import json
        from scipy.interpolate import splev

        coef_v = json.loads(self.coef_v)
        coef_u = json.loads(self.coef_u)

        v = splev(T, coef_v)
        u = splev(T, coef_u)
        return float(v), float(u)
