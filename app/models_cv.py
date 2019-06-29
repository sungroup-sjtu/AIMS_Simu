import json
from sqlalchemy import Integer, Text, Column, DateTime
from datetime import datetime
from . import db
from .models import NotNullColumn


class Cv(db.Model):
    __bind_key__ = 'cv'
    __tablename__ = 'cv'
    id = NotNullColumn(Integer, primary_key=True)
    name = NotNullColumn(Text)
    smiles = NotNullColumn(Text)
    post_result = NotNullColumn(Text)
    time = Column(DateTime, default=datetime.now)

    def __repr__(self):
        return '<Cv: %i: %s>' % (self.id, self.smiles)

    def get_post_cv(self, T=298):
        from numpy.polynomial.polynomial import polyval
        coef = json.loads(self.post_result).get('Cv')
        return polyval(T, coef)

    def get_post_enthalpy(self, T=298):
        from numpy.polynomial.polynomial import polyval
        coef = json.loads(self.post_result).get('enthalpy')
        return polyval(T, coef) * 627.51 * 4.184

    @staticmethod
    def load_from_log(log_cv, log_enthalpy):
        with open(log_cv) as f:
            lines = f.read().splitlines()
        for line in lines:
            if line == '' or line.startswith('#'):
                continue
            if 'failed' in line or 'Err' in line:
                continue

            words = line.strip().split()
            name = words[1]
            smiles = words[2].replace('>', '')
            coef = list(map(float, words[3:8]))
            score = float(words[8])

            def get_atom_number(smiles):
                import pybel
                py_mol = pybel.readstring('smi', smiles)
                py_mol.addh()
                return len(py_mol.atoms)

            if score < 0.999 and not(get_atom_number(smiles) == 1 and score == 0.0):
                print(name, smiles, 'Score < 0.999')
                continue

            count = Cv.query.filter(Cv.smiles == smiles).count()
            if count == 0:
                dict = {'Cv': coef}
                cv = Cv(name=name, smiles=smiles, post_result=json.dumps(dict))
                db.session.add(cv)
            elif count == 1:
                cv = Cv.query.filter(Cv.smiles == smiles).first()
                result = json.loads(cv.post_result)
                if result.get('Cv') == None:
                    dict = {'Cv': coef}
                    result.update(dict)
                cv.post_result = json.dumps(result)
            else:
                raise Exception('Molecule %s exists %i times in Cv database' % (smiles, count))
        db.session.commit()

        with open(log_enthalpy) as f:
            lines = f.read().splitlines()
        for line in lines:
            if line == '' or line.startswith('#'):
                continue
            if 'failed' in line or 'Err' in line:
                continue

            words = line.strip().split()
            name = words[1]
            smiles = words[2].replace('>', '')
            coef = list(map(float, words[3:6]))
            score = float(words[6])

            if score < 0.999:
                print(name, smiles, 'Score < 0.999')
                continue

            count = Cv.query.filter(Cv.smiles==smiles).count()
            if count==0:
                dict = {'enthalpy': coef}
                cv = Cv(name=name, smiles=smiles, post_result=json.dumps(dict))
                db.session.add(cv)
            elif count==1:
                cv = Cv.query.filter(Cv.smiles==smiles).first()
                result = json.loads(cv.post_result)
                if result.get('enthalpy') == None:
                    dict = {'enthalpy': coef}
                    result.update(dict)
                cv.post_result = json.dumps(result)
            else:
                raise Exception('Molecule %s exists %i times in Cv database' % (smiles, count))
        db.session.commit()
