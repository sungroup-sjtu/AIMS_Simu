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

    def get_post_result(self, T=298):
        from numpy.polynomial.polynomial import polyval
        coef = json.loads(self.post_result)
        return polyval(T, coef)
