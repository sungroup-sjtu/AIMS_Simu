from flask import request, jsonify
from . import api


@api.route('/api')
def index():
    data = {'ip': request.remote_addr}
    return jsonify(data)


@api.route('/api/check', methods=['GET', 'POST'])
def check():
    data = {'success': True, 'ip': request.remote_addr}
    return jsonify(data)


@api.route('/api/submit', methods=['GET', 'POST'])
def submit():
    data = {'success': False, 'ip': request.remote_addr}
    return jsonify(data)
