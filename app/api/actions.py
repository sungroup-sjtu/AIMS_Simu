import json

from flask import request

from mstools.utils import get_T_list_from_range, get_P_list_from_range
from ..models import *


class ComputeAction():
    def init_from_json(self, compute_json) -> int:
        detail_json = compute_json['detail']
        n_components = len(detail_json['smiles_list'])
        for smiless in detail_json['smiles_list']:
            if len(smiless) == 0:
                raise Exception('Invalid information')
        procedures = detail_json['procedures']
        if len(procedures) == 0:
            raise Exception('Invalid information')

        compute = Compute()
        compute.web_id = compute_json['id']
        compute.web_user_id = compute_json['user_id']
        compute.web_ip = request.remote_addr
        compute.n_components = n_components
        compute.json = json.dumps(detail_json)

        db.session.add(compute)
        db.session.flush()

        T_list = get_T_list_from_range(detail_json['t_min'], detail_json['t_max'])
        P_list = get_P_list_from_range(detail_json['p_min'], detail_json['p_max'])

        if n_components == 1:
            combinations = []
            for smiles in detail_json['smiles_list'][0]:
                for procedure in procedures:
                    for T in T_list:
                        if procedure not in ComputeProcedure.T_RELEVANT:
                            T = None
                        for P in P_list:
                            if procedure not in ComputeProcedure.P_RELEVANT:
                                P = None

                            combination = [smiles, procedure, T, P]
                            if combination in combinations:
                                continue

                            combinations.append(combination)
                            job = Task()
                            job.compute_id = compute.id
                            job.smiles = smiles
                            job.procedure = procedure
                            job.t = T
                            job.p = P
                            db.session.add(job)

        elif n_components == 2:
            for smiles in detail_json['smiles_list'][0]:
                for smiles2 in detail_json['smiles_list'][1]:
                    for procedure in procedures:
                        for T in T_list:
                            for P in P_list:
                                job = JobBinary()
                                job.compute_id = compute.id
                                job.smiles = smiles
                                job.smiles2 = smiles2
                                job.procedure = procedure
                                job.t = T
                                job.p = P
                                db.session.add(job)

        try:
            db.session.commit()
        except Exception as e:
            db.session.rollback()
            raise Exception(str(e))

        job_thread = JobThread(compute.id)
        job_thread.start()

        return compute.id
