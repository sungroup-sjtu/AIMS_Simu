#!/usr/bin/env python3
# coding=utf-8

import sys
sys.path.append('../AIMS_Simu')
from app import create_app
from app.models import *
from app.models_nist import *

procedure = 'npt'
app = create_app(procedure)
app.app_context().push()

datas = NistData.query
sp1 = NistSpline.query.filter(NistSpline.molecule_id == 4).filter(NistSpline.property_id == 1).first()
sp2 = NistSpline.query.filter(NistSpline.molecule_id == 4).filter(NistSpline.property_id == 3).first()