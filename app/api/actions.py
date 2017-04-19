from flask import request

from ..models import *


class ComputeAction():
    def init_from_json(self, compute_json) -> int:
        detail_json = compute_json['detail']

        compute = Compute()
        compute.web_id = compute_json['id']
        compute.web_user_id = compute_json['user_id']
        compute.web_ip = request.remote_addr
        compute.json = json.dumps(detail_json)

        db.session.add(compute)
        db.session.flush()

        try:
            compute.create_tasks()
            db.session.commit()

        except Exception as e:
            db.session.rollback()
            raise Exception(str(e))

        # job_thread = TaskThread(current_app._get_current_object(), compute.id)
        # job_thread.start()

        return compute.id
