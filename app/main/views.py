from flask import render_template, abort
from . import main
from ..models import Compute, Task, Job


@main.route('/')
def index():
    return render_template('index.html')


@main.route('/compute/<int:compute_id>')
def show_compute(compute_id):
    compute = Compute.query.get(compute_id)
    return render_template('compute.html', compute=compute)


@main.route('/task/<int:task_id>')
def show_task(task_id):
    task = Task.query.get(task_id)
    return render_template('task.html', task=task)


@main.route('/job/<int:job_id>')
def show_job(job_id):
    job = Job.query.get(job_id)
    return render_template('job.html', job=job)
