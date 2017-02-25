from flask import request

from . import api
from .actions import *

import json


@api.route('/')
def index():
    data = {'ip': request.remote_addr}
    return json.dumps(data)


@api.route('/check', methods=['POST'])
def check():
    compute_id = request.form['compute_id']
    user_id = request.form['user_id']
    compute = Compute.query.filter(and_(Compute.id == compute_id, Compute.web_user_id == user_id)).first()
    if compute == None:
        return json.dumps({'success': False,
                           'reason': 'Not found'})

    if compute.type == Compute.Type.UNARY:
        tasks = ComputeUnary.query.filter(ComputeUnary.compute_id == compute.id).all()
    elif compute.type == Compute.Type.BINARY:
        tasks = ComputeBinary.query.filter(ComputeBinary.compute_id == compute.id).all()

    data = {'success': True, 'ip': request.remote_addr}
    return json.dumps(data)


@api.route('/submit', methods=['POST'])
def submit():
    compute_json = request.json

    print(compute_json)

    computeAction = ComputeAction()
    try:
        compute_id = computeAction.init_from_json(compute_json)
    except Exception as e:
        return json.dumps({'success': False,
                           'ip': request.remote_addr,
                           'reason': str(e)
                           })
    else:
        return json.dumps({'success': True,
                           'compute_id': compute_id,
                           'ip': request.remote_addr
                           })
