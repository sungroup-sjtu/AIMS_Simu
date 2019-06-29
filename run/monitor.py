#!/usr/bin/env python3
# coding=utf-8

import os, sys, time
from sqlalchemy import func, or_

sys.path.append('..')
from app import create_app
from app.models import *
# necessary functions
def detect_exit():
    os.chdir(CWD)
    if os.path.exists('EXIT-' + procedure):
        print('EXIT file detected')
        os.remove('EXIT-' + procedure)
        sys.exit()
def process_pbs_job(n_pbs=20):
    for pbs_job in PbsJob.query.filter(PbsJob.submitted == False).limit(n_pbs):
        detect_exit()
        pbs_job.submit()
def process_task_build(n_task=20, random=False):
    tasks = Task.query.filter(Task.stage == Compute.Stage.SUBMITTED).filter(Task.status == Compute.Status.DONE).filter(
        Task.procedure == procedure)
    if random:
        tasks = tasks.order_by(func.random())
    for task in tasks.limit(n_task):
        detect_exit()
        task.build()
        if task.stage == Compute.Stage.BUILDING and task.status == Compute.Status.DONE:
            task.run()
def process_task_run(n_task=20, random=False):
    tasks = Task.query.filter(Task.stage == Compute.Stage.BUILDING).filter(Task.status == Compute.Status.DONE).filter(
        Task.procedure == procedure)
    if random:
        tasks = tasks.order_by(func.random())
    for task in tasks.limit(n_task):
        detect_exit()
        if task.run() == -1:
            break
def process_task_check(n_task=20):
    app.jobmanager.update_stored_jobs()
    app.jm_extend.update_stored_jobs()
    tasks = Task.query.filter(Task.stage == Compute.Stage.RUNNING).filter(Task.status == Compute.Status.STARTED) \
        .filter(Task.procedure == procedure).limit(n_task)
    for task in tasks:
        detect_exit()
        task.check_finished_multiprocessing()
def process_extend(n_job=None):
    jobs_extend_tmp = Job.query.filter(Job.status == Compute.Status.ANALYZED).filter(Job.converged == False) \
        .filter(Job.cycle < Config.EXTEND_CYCLE_LIMIT)
    jobs_extend = []
    for job in jobs_extend_tmp:
        if job.task.procedure == procedure and job.result is not None:
            jobs_extend.append(job)
    if n_job is not None and len(jobs_extend) > n_job:
        jobs_extend = jobs_extend[0:n_job]
    current_app.logger.info('Extend all jobs')
    if not current_app.jm_extend.is_working():
        current_app.logger.warning('JobManager not working')
        return

    n_pbs_extend = len(jobs_extend)
    if current_app.config['EXTEND_GMX_MULTI']:
        # do not extend when few jobs need extend
        if len(jobs_extend) < current_app.config['EXTEND_GMX_MULTI_NJOB']:
            current_app.logger.warning('%i job need extend' % (len(jobs_extend)))
            return

        name_list = json.loads(jobs_extend[0].result).get('name')
        extend_jobs_dict = {}
        for name in name_list:
            extend_jobs_dict.update({name: {}})
        for job in jobs_extend:
            result = json.loads(job.result)
            continue_list = result.get('continue')
            continue_n = result.get('continue_n')
            for i, name in enumerate(name_list):
                if continue_list[i] == True:
                    simulation_part_n = job.get_simulation_part(name)
                    if extend_jobs_dict.get(name).get(str(simulation_part_n)) == None:
                        extend_jobs_dict.get(name)[str(simulation_part_n)] = []
                    extend_jobs_dict.get(name).get(str(simulation_part_n)).append(ExtendJob(job, continue_n[i]))
        n_pbs_extend = 0
        multi_njob = current_app.config['EXTEND_GMX_MULTI_NJOB']
        for name in name_list:
            for simulation_part_n in extend_jobs_dict.get(name).keys():
                extend_jobs_dict.get(name).get(simulation_part_n).sort(reverse=True)
                n_pbs = math.floor(len(extend_jobs_dict.get(name).get(simulation_part_n)) / multi_njob)
                n_jobs = n_pbs * multi_njob
                extend_jobs_dict.get(name)[simulation_part_n] = extend_jobs_dict.get(name).get(simulation_part_n)[
                                                                0:n_jobs]
                n_pbs_extend += n_pbs
                for i in range(n_pbs):
                    continue_n = extend_jobs_dict.get(name).get(simulation_part_n)[i * multi_njob].continue_n
                    for j in range(1, multi_njob):
                        n_job = i * multi_njob + j
                        extend_jobs_dict.get(name).get(simulation_part_n)[n_job].continue_n = continue_n
                current_app.logger.info('%s-%i: %i jobs need extend' % (name, int(simulation_part_n), n_jobs))

    if current_app.jm_extend.n_running_jobs > current_app.config['EXTEND_PBS_NJOB_LIMIT']:
        current_app.logger.warning('EXTEND_PBS_NJOB_LIMIT reached')
        return

    try:
        if current_app.config['EXTEND_GMX_MULTI']:
            job_list = []
            sim = init_simulation(procedure, extend=True)
            for name in name_list:
                for simulation_part_n in extend_jobs_dict.get(name).keys():
                    multi_dirs = []
                    if extend_jobs_dict.get(name).get(simulation_part_n) == []:
                        continue
                    for a in extend_jobs_dict.get(name).get(simulation_part_n):
                        job = a.job
                        if job not in job_list:
                            job.cycle += 1
                            job.status = Compute.Status.STARTED
                            job.pbs_jobs_id = json.dumps([])
                            job_list.append(job)
                        continue_n = a.continue_n
                        os.chdir(job.dir)
                        multi_dirs.append(job.dir)
                        multi_cmds = sim.extend_single(jobname='%s-%s-%i' % (job.name, name, job.cycle + 1),
                                                       sh='_job.extend-%s-%i.sh' % (name, job.cycle + 1), name=name,
                                                       continue_n=continue_n)

                    commands_list = GMX.generate_gpu_multidir_cmds(multi_dirs, multi_cmds,
                                                                   n_parallel=multi_njob,
                                                                   n_gpu=current_app.jm_extend.ngpu,
                                                                   n_procs=current_app.jm_extend.nprocs)
                    extend_dir = os.path.join(current_app.config['WORK_DIR'], procedure, 'extend_sh')
                    if not os.path.exists(extend_dir):
                        os.mkdir(extend_dir)
                    os.chdir(extend_dir)

                    n = 1
                    while True:
                        sh = os.path.join(extend_dir, '_job.extend-%s-%i.sh' % (name, n))
                        if os.path.exists(sh):
                            n += 1
                        else:
                            break

                    for i, commands in enumerate(commands_list):
                        sh = os.path.join(extend_dir, '_job.extend-%s-%i.sh' % (name, n + i))
                        if procedure == 'npt-2':
                            pbs_name = '%s-2-global-extend-%i' % (name, n + i)
                        elif procedure == 'npt-3':
                            pbs_name = '%s-3-global-extend-%i' % (name, n + i)
                        else:
                            pbs_name = '%s-global-extend-%i' % (name, n + i)

                        current_app.jm_extend.generate_sh(extend_dir, commands, name=pbs_name, sh=sh)

                        # instead of run directly, we add a record to pbs_job
                        pbs_job = PbsJob(extend=True)
                        pbs_job.name = pbs_name
                        pbs_job.sh_file = sh
                        db.session.add(pbs_job)
                        db.session.flush()

                        for a in extend_jobs_dict.get(name).get(simulation_part_n)[i * multi_njob:(i + 1) * multi_njob]:
                            job = a.job
                            pbs_jobs_id = json.loads(job.pbs_jobs_id)
                            pbs_jobs_id.append(pbs_job.id)
                            job.pbs_jobs_id = json.dumps(pbs_jobs_id)
                        db.session.commit()

                        # submit job, record if success or failed
                        pbs_job.submit()
                        time.sleep(0.2)
        else:
            for job in jobs_extend:
                job.extend()
                time.sleep(0.2)



    except Exception as e:
        current_app.logger.error('Extend jobs failed %s ' % (repr(e)))
        traceback.print_exc()
        return 0
    else:
        return n_pbs_extend
def process_post_process():
    current_app.logger.info('Post process finished tasks')
    tasks = Task.query.filter(Task.procedure == procedure).filter(Task.post_result == None)
    current_app.logger.info('Total task = %i' % (tasks.count()))
    count = 1
    for task in tasks:
        sys.stdout.write('\rpresent task = %i' % (count))
        count += 1
        task.post_process(overwrite=False)
        task.get_LJ_atom_type()

def config_check():
    if Config.LJ96:
        if current_app.config['PBS_ARGS'][0] == 'gtx' or current_app.config['EXTEND_PBS_ARGS'][0] == 'gtx':
            print('GPU cannot use for LJ96 job')
            return False
    return True
# in ppm simulation, in a task, if the viscosity of highest temperature
def extend_unfinished_jobs():
    current_app.logger.info('Continue abnormally terminated ppm jobs')
    tasks = Task.query.filter(Task.procedure == procedure)
    for task in tasks:
        for job in task.jobs:
            job.continue_terminated_job()
def process_task_bug_fix():
    # In npt simulation. sometimes the simulation will ended abnormally due to unstable initial structure,
    # this problem may be fixed by rebuild the system using large system
    '''
    tasks = Task.query.filter(Task.status == Compute.Status.FAILED)
    for task in tasks:
        job_status = []
        job_result = []
        for job in task.jobs:
            job_status.append(job.status)
            job_result.append(job.result)
        if set(job_status) == {Compute.Status.FAILED} and set(job_result) == {None}:
            if task.procedure == 'npt':
                task.reset(add_mol_list=[20] * task.n_components)  # add 20 molecules of each component
            elif task.procedure == 'ppm':
                task.reset()

    if procedure == 'ppm':
        jobs = Job.query.filter(Job.status == Compute.Status.FAILED).filter(Job.result == None)
        for job in jobs:
            if job.task.procedure == 'ppm':
                job.status = Compute.Status.STARTED
                if job.task.status == Compute.Status.DONE:
                    job.task.status = Compute.Status.STARTED
    db.session.commit()
    # there is a bug when analyze the simulation result in ppm calculation,
    # '(sqlite3.DatabaseError) database disk image is malformed'
    # so far cannot fix it directly, using following function to reanalyze the result
    '''
# this function is used, when you increase the REPEAT_NUMBER in config.py, to generate more repeated jobs with
# different random seed
def check_tasks_jobs():
    tasks = Task.query.filter(Task.procedure == procedure)
    for task in tasks:
        task.check_task_jobs()


import argparse
parser = argparse.ArgumentParser(description='This is a code to monitor the high-throughput simulation')
parser.add_argument('-p', '--procedure', type=str, help='procedure of the compute: npt(ppm), or nvt-slab')
opt = parser.parse_args()

procedure = opt.procedure

app = create_app(procedure)
app.app_context().push()

CWD = os.getcwd()
while True:
    # process_task_bug_fix()
    if not config_check():
        break
    # check_tasks_jobs() # use when you increase the repeat number in config.py

    process_pbs_job(n_pbs=50)
    process_task_run(n_task=20)
    process_task_build(n_task=20)
    process_task_check(n_task=4000)
    if procedure == 'ppm':
        extend_unfinished_jobs()  # in ppm simulation, in a task, if the viscosity of highest temperature is too slow, then the
    # simulation failed in strong acceleration condition. Then the simulation at low temperature will also failed.
    # When GPU is used. single simulation failure will terminate the whole GPU task, use extend_unfinished_jobs() to fix this problem.
    process_extend(n_job=200)
    process_post_process()

    app.logger.info('Sleep 1800 seconds ...')
    time.sleep(1800)

