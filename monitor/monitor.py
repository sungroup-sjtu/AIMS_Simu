from app.models import *

tasks = Task.query.all()
for task in tasks:
    if task.stage == Compute.Stage.SUBMITTED and task.status == Compute.Status.DONE:
        try:
            task.build()
        except Exception as e:
            print('Error when build task %s: %s' % (repr(task), str(e)))

    elif task.stage == Compute.Stage.BUILDING and task.status == Compute.Status.DONE:
        try:
            task.run()
        except Exception as e:
            print('Error when run task %s: %s' % (repr(task), str(e)))

    elif task.stage == Compute.Stage.RUNNING and task.stage == Compute.Status.STARTED:
        for job in task.jobs:
            job.check()
