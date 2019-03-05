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

    def get_post_data(self, T=298):
        from numpy.polynomial.polynomial import polyval
        coef = json.loads(self.post_result)
        return polyval(T, coef)

    @staticmethod
    def load_from_log(log):
        with open(log) as f:
            lines = f.read().splitlines()

        cv_list = []
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

            if score < 0.999:
                print(name, smiles, 'Score < 0.999')
                continue

            cv = Cv(name=name, smiles=smiles, post_result=json.dumps(coef))
            cv_list.append(cv)

        print('%i rows inserted' % len(cv_list))
        db.session.add_all(cv_list)
        db.session.commit()
