import pybel
from sqlalchemy import Integer, Text, Column, Boolean, String, UniqueConstraint, ForeignKey, Float
from . import db
from sqlalchemy.orm import relationship
from scipy import interpolate
import json


class Property(db.Model):
    __bind_key__ = 'ilthermo'
    __tablename__ = 'property'
    id = Column(Integer, primary_key=True)
    name = Column(Text, unique=True)
    unit = Column(Text)

    datas = relationship('Data', lazy='dynamic')

    def __repr__(self):
        return '<Property: %i %s>' % (self.id, self.name)


class Paper(db.Model):
    __bind_key__ = 'ilthermo'
    __tablename__ = 'paper'
    id = Column(Integer, primary_key=True)
    year = Column(Integer)
    title = Column(Text, unique=True)
    author = Column(Text)

    datas = relationship('Data', lazy='dynamic')

    def __repr__(self):
        return '<Paper: %i %s>' % (self.id, self.author)


class DataSet(db.Model):
    __bind_key__ = 'ilthermo'
    __tablename__ = 'dataset'
    id = Column(Integer, primary_key=True)
    code = Column(String(5))
    searched = Column(Boolean)


class Ion(db.Model):
    __bind_key__ = 'ilthermo'
    __tablename__ = 'ion'
    id = Column(Integer, primary_key=True)
    charge = Column(Integer)
    name = Column(Text, unique=True)
    searched = Column(Boolean)
    popular = Column(Boolean, default=False)
    selected = Column(Boolean, default=False)
    smiles = Column(Text)
    iupac = Column(Text)
    ignored = Column('validated', Boolean, default=False)
    category = Column(Text)
    n_paper = Column(Integer)
    times = Column(Integer)
    category_xy = Column(Integer)
    force_field_support = Column(Integer)

    molecules_cation = relationship('Molecule', lazy='dynamic', foreign_keys='Molecule.cation_id')
    molecules_anion = relationship('Molecule', lazy='dynamic', foreign_keys='Molecule.anion_id')

    def __repr__(self):
        if self.charge > 0:
            return '<Ion +%i: %i %s>' % (abs(self.charge), self.id, self.name)
        else:
            return '<Ion -%i: %i %s>' % (abs(self.charge), self.id, self.name)

    @property
    def molecules(self):
        if self.charge > 0:
            return self.molecules_cation
        else:
            return self.molecules_anion

    @property
    def n_heavy(self):
        try:
            py_mol = pybel.readstring('smi', self.smiles)
        except:
            raise Exception('Smiles not valid')
        return py_mol.OBMol.NumHvyAtoms()

    def update_formula(self):
        self.formula = pybel.readstring('smi', self.smiles).formula


class Molecule(db.Model):
    __bind_key__ = 'ilthermo'
    __tablename__ = 'molecule'
    __table_args__ = (UniqueConstraint('cation_id', 'anion_id', name='ion_id'),)
    id = Column(Integer, primary_key=True)
    code = Column(String(6))
    name = Column(Text, unique=True)
    cation_id = Column(Integer, ForeignKey(Ion.id))
    anion_id = Column(Integer, ForeignKey(Ion.id))
    formula = Column(Text)
    popular = Column(Boolean, default=False)
    selected = Column(Boolean, default=False)
    fit = Column(Text)
    cation_category = Column(Integer)
    anion_category = Column(Integer)

    cation = relationship('Ion', foreign_keys='Molecule.cation_id')
    anion = relationship('Ion', foreign_keys='Molecule.anion_id')

    datas = relationship('Data', lazy='dynamic')
    splines = relationship('Spline', lazy='dynamic')

    def __repr__(self):
        return '<Molecule: %i %s>' % (self.id, self.name)

    def smiles(self):
        if self.cation.smiles is None or self.anion.smiles is None:
            return None  # '%i.%i' % (self.cation.id, self.anion.id)
        else:
            return self.cation.smiles + '.' + self.anion.smiles


class Data(db.Model):
    __bind_key__ = 'ilthermo'
    __tablename__ = 'data'
    id = Column(Integer, primary_key=True)
    molecule_id = Column(Integer, ForeignKey(Molecule.id))
    paper_id = Column(Integer, ForeignKey(Paper.id))
    property_id = Column(Integer, ForeignKey(Property.id))
    phase = Column(String(20))
    t = Column(Float)
    p = Column(Float, nullable=True)
    value = Column(Float)
    stderr = Column(Float)

    molecule = relationship('Molecule', foreign_keys='Data.molecule_id')
    paper = relationship('Paper', foreign_keys='Data.paper_id')
    property = relationship('Property', foreign_keys='Data.property_id')

    def __repr__(self):
        return '<Data: %s: %.1f %.1f %f>' % (self.property.name, self.t or 0, self.p or 0, self.value)


class Spline(db.Model):
    __bind_key__ = 'ilthermo'
    __tablename__ = 'spline'
    id = Column(Integer, primary_key=True)
    molecule_id = Column(Integer, ForeignKey(Molecule.id))
    property_id = Column(Integer, ForeignKey(Property.id))
    paper_id = Column(Integer, ForeignKey(Paper.id))
    t_min = Column(Float)
    t_max = Column(Float)
    coef_v = Column(Text)  # value
    coef_u = Column(Text)  # uncertainty
    coef_VTF = Column(Text)

    molecule = relationship(Molecule)
    property = relationship(Property)
    paper = relationship(Paper)

    def get_data(self, T):
        if T < self.t_min or T > self.t_max:
            return None, None

        coef_v = json.loads(self.coef_v)
        coef_u = json.loads(self.coef_u)

        v = interpolate.splev(T, coef_v)
        u = interpolate.splev(T, coef_u)
        return float(v), float(u)

    def get_VTF_data(self, T):
        if T < self.t_min or T > self.t_max:
            return None

        from mstools.analyzer.fitting import VTFval
        coef_VTF, score = json.loads(self.coef_VTF)
        if score > 0.95:
            return VTFval(T, coef_VTF)
        else:
            return None
