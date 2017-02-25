from . import db
from sqlalchemy import Column, ForeignKey, Integer, Text, String, Boolean, DateTime
from datetime import datetime

from functools import partial

NotNullColumn = partial(Column, nullable=False)


class Simulation(db.Model):
    __tablename__ = 'simulation'
    id = NotNullColumn(Integer, primary_key=True)
    name = NotNullColumn(String(200))

    properties = db.relationship('Property', lazy='dynamic')

    class Type:
        NPT = 1
        NVT = 2
        SLAB = 3
        VISCOSITY = 4
        BINARY_SLAB = 5
        SOLVATION_FE = 6

    @staticmethod
    def insert_simulations():
        simulations = [
            (Simulation.Type.NPT, 'npt'),
            (Simulation.Type.NVT, 'nvt'),
            (Simulation.Type.SLAB, 'slab'),
            (Simulation.Type.VISCOSITY, 'viscosity'),
            (Simulation.Type.BINARY_SLAB, 'binary_slab'),
            (Simulation.Type.SOLVATION_FE, 'solvation_fe'),
        ]
        for i in simulations:
            simulation = Simulation(id=i[0], name=i[1])
            db.session.add(simulation)
        try:
            db.session.commit()
        except Exception as e:
            print(str(e))
            db.session.rollback()


class Property(db.Model):
    __tablename__ = 'property'
    id = NotNullColumn(Integer, primary_key=True)
    name = NotNullColumn(db.String(200), unique=True)
    abbr = NotNullColumn(db.String(200), unique=True)
    unit = NotNullColumn(db.String(200))
    type = NotNullColumn(db.Integer)
    simulation_id = NotNullColumn(Integer, ForeignKey(Simulation.id))

    simulation = db.relationship(Simulation)

    def __repr__(self):
        return '<Property: %s>' % self.name

    class Type:
        TP = 1
        SATURATION = 2
        CRITICAL = 3
        choices = [TP, SATURATION, CRITICAL]
        T_RELEVANT = [TP, SATURATION]
        P_RELEVANT = [TP]

    @staticmethod
    def insert_properties():
        properties = [
            # TP thermodynamic properties
            ('Density', 'density', 'g/mL', Property.Type.TP, Simulation.Type.NPT),
            ('Enthalpy of vaporization', 'HoV', 'kcal/mol', Property.Type.TP, Simulation.Type.NPT),
            ('Isochoric heat capacity', 'Cv', 'cal/mol/K', Property.Type.TP, Simulation.Type.NPT),
            ('Isobaric heat capacity', 'Cp', 'cal/mol/K', Property.Type.TP, Simulation.Type.NPT),
            ('Isothermal compressibility', 'compressibility', '/Pa', Property.Type.TP, Simulation.Type.NPT),
            ('Volume expansivity', 'expansivity', '/K', Property.Type.TP, Simulation.Type.NPT),
            ('Joule-Thomson coefficient', 'JT', 'K/Pa', Property.Type.TP, Simulation.Type.NPT),
            ('Speed of sound', 'cSound', 'm/s', Property.Type.TP, Simulation.Type.NPT),

            # TP dynamic properties
            ('Diffusion coefficient', 'diffusion', 'm^2/s', Property.Type.TP, Simulation.Type.NVT),
            ('Viscosity', 'viscosity', 'cP', Property.Type.TP, Simulation.Type.VISCOSITY),
            ('Thermal conductivity', 'TC', 'W/m/K', Property.Type.TP, Simulation.Type.NVT),

            # Saturation properties
            ('Surface tension', 'ST', 'mN/m', Property.Type.SATURATION, Simulation.Type.SLAB),
            ('Liquid density', 'dLiquid', 'g/mL', Property.Type.SATURATION, Simulation.Type.SLAB),
            ('Vapor density', 'dVapor', 'g/mL', Property.Type.SATURATION, Simulation.Type.SLAB),
            ('Vapor pressure', 'pVapor', 'kPa', Property.Type.SATURATION, Simulation.Type.SLAB),

            # Critical properties
            ('Critical temperature', 'tCritical', 'K', Property.Type.CRITICAL, Simulation.Type.SLAB),
            ('Critical pressure', 'pCritical', 'kPa', Property.Type.CRITICAL, Simulation.Type.SLAB),
            ('Critical density', 'dCritical', 'g/mL', Property.Type.CRITICAL, Simulation.Type.SLAB),

            # Binary TP properties
            ('Solvation free energy', 'gSolvation', 'kcal/mol', Property.Type.TP, Simulation.Type.SOLVATION_FE),
            ('Solubility', 'solubility', 'mol/L', Property.Type.TP, Simulation.Type.BINARY_SLAB),

            # Binary Saturation properties
            ('Mole fraction in gas phase', 'Xg', '', Property.Type.SATURATION, Simulation.Type.BINARY_SLAB),
            ('Mole fraction in liquid phase', 'Xl', '', Property.Type.SATURATION, Simulation.Type.BINARY_SLAB),
        ]

        for i in properties:
            property = Property(name=i[0], abbr=i[1], unit=i[2], type=i[3], simulation_id=i[4])
            db.session.add(property)
        try:
            db.session.commit()
        except Exception as e:
            print(str(e))
            db.session.rollback()


class Compute(db.Model):
    __tablename__ = 'compute'
    id = NotNullColumn(Integer, primary_key=True)
    web_id = NotNullColumn(Integer)
    web_user_id = NotNullColumn(Integer)
    web_ip = NotNullColumn(String(200))
    type = NotNullColumn(Integer)
    time = NotNullColumn(DateTime, default=datetime.now)

    class Type:
        UNARY = 1
        BINARY = 2

    class Status:
        SUBMITTED = 0
        PREPARING = 1
        RUNNING = 2
        DONE = 9
        FAILED = -1

        text = {
            SUBMITTED: 'Submitted',
            PREPARING: 'Preparing',
            RUNNING: 'Running',
            DONE: 'Done',
            FAILED: 'Failed'
        }


class ComputeUnary(db.Model):
    __tablename__ = 'compute_unary'
    id = NotNullColumn(Integer, primary_key=True)
    compute_id = NotNullColumn(Integer, ForeignKey(Compute.id))
    web_molecule_id = NotNullColumn(Integer)
    web_molecule_smiles = NotNullColumn(Text)
    simulation_id = NotNullColumn(Integer, ForeignKey(Simulation.id))
    t_min = NotNullColumn(Integer)
    t_max = NotNullColumn(Integer)
    p_min = NotNullColumn(Integer)
    p_max = NotNullColumn(Integer)
    time = NotNullColumn(DateTime, default=datetime.now)
    job_name = Column(String(200), nullable=True)
    status = NotNullColumn(Integer, default=Compute.Status.SUBMITTED)
    remark = Column(Text, nullable=True)

    compute = db.relationship(Compute)
    simulation = db.relationship(Simulation)

    def generate_tasks(self):
        pass

    def build(self):
        pass

    def run(self):
        pass

    def get_dir(self):
        pass

    def analyze(self):
        pass


class ComputeBinary(db.Model):
    __tablename__ = 'compute_binary'
    id = NotNullColumn(Integer, primary_key=True)
    compute_id = NotNullColumn(Integer, ForeignKey(Compute.id))
    web_id = NotNullColumn(Integer)
    web_molecule_id = NotNullColumn(Integer)
    web_molecule_smiles = NotNullColumn(Text)
    web_molecule2_id = NotNullColumn(Integer)
    web_molecule2_smiles = NotNullColumn(Text)
    simulation_id = NotNullColumn(Integer, ForeignKey(Simulation.id))
    t_min = NotNullColumn(Integer)
    t_max = NotNullColumn(Integer)
    p_min = NotNullColumn(Integer)
    p_max = NotNullColumn(Integer)
    time = NotNullColumn(DateTime, default=datetime.now)
    name = Column(String(200), nullable=True)
    status = NotNullColumn(Integer, default=Compute.Status.SUBMITTED)
    remark = Column(Text, nullable=True)

    compute = db.relationship(Compute)
    simulation = db.relationship(Simulation)
