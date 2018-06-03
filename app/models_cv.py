import json
from sqlalchemy import Integer, Text
from . import db
from .models import NotNullColumn


class Cv(db.Model):
    __bind_key__ = 'cv'
    __tablename__ = 'cv'
    id = NotNullColumn(Integer, primary_key=True)
    cas = NotNullColumn(Text)
    smiles = NotNullColumn(Text)
    post_result = NotNullColumn(Text)

    def __repr__(self):
        return '<Cv: %i: %s %s>' % (self.id, self.cas, self.smiles)

    def get_post_data(self, T=298):
        from numpy.polynomial.polynomial import polyval
        coef = json.loads(self.post_result)
        return polyval(T, coef)

    @staticmethod
    def load_from_log(log):
        with open(log) as f:
            lines = f.read().splitlines()

        for line in lines:
            if line == '' or line.startswith('#'):
                continue
            if 'failed' in line or 'Errno' in line:
                continue

            words = line.strip().split()
            cas = words[1]
            smiles = words[2]
            smiles = smiles.replace('>', '')
            coef = list(map(float, words[3:8]))
            cv = Cv(cas=cas, smiles=smiles, post_result=json.dumps(coef))
            db.session.add(cv)
        try:
            db.session.commit()
        except Exception as e:
            print(repr(e))
            db.session.rollback()
