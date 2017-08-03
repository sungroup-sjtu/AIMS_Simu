from flask import request

from ..models import *


class ComputeAction():
    def init_from_json_dict(self, json_dict) -> int:
        compute = Compute()
        compute.web_ip = request.remote_addr
        try:
            compute.web_id = json_dict['id']
            compute.web_user_id = json_dict['user_id']
            compute.json = json.dumps(json_dict['detail'])
            compute.remark = json_dict['remark']
        except:
            raise Exception('Incomplete information')

        db.session.add(compute)
        db.session.commit()

        try:
            compute.create_tasks()
        except Exception:
            db.session.delete(compute)
            db.session.commit()
            raise

        return compute.id
