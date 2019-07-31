#!/usr/bin/env python3
# coding=utf-8

import os, sys, time

sys.path.append('..')
from app import create_app
from app.models import Task, Compute, PbsJob

app = create_app(sys.argv[1])
app.app_context().push()


def main():
    tasks = Task.query
    n_total = tasks.count()
    for i, task in enumerate(tasks):
        print(f'{i} / {n_total}', task.post_result(force=True))


if __name__ == '__main__':
    main()
