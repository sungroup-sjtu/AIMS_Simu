import sys
import re
import json
import shutil
import time
from datetime import datetime
from functools import partial

from sqlalchemy import Column, ForeignKey, Integer, Text, String, Boolean, DateTime, and_

from . import db, log

from config import Config

sys.path.append(Config.MS_TOOLS_DIR)
from mstools.simulation.procedure import Procedure
from mstools.utils import *

NotNullColumn = partial(Column, nullable=False)

from mstools.jobmanager import *
from mstools.wrapper import GMX

if Config.PBS_MANAGER == 'local':
    jobmanager = Local(queue_dict=Config.PBS_QUEUE_DICT, env_cmd=Config.PBS_ENV_CMD)
elif Config.PBS_MANAGER == 'torque':
    jobmanager = Torque(queue_dict=Config.PBS_QUEUE_DICT, env_cmd=Config.PBS_ENV_CMD)
elif Config.PBS_MANAGER == 'slurm':
    jobmanager = Slurm(queue_dict=Config.PBS_QUEUE_DICT, env_cmd=Config.PBS_ENV_CMD)
else:
    raise Exception('Job manager not supported')


def init_simulation(procedure):
    from mstools.simulation import gmx as simulationEngine
    kwargs = {'packmol_bin': Config.PACKMOL_BIN, 'dff_root': Config.DFF_ROOT,
              'gmx_bin': Config.GMX_BIN, 'jobmanager': jobmanager}

    if procedure == 'npt':
        return simulationEngine.Npt(**kwargs)
    elif procedure == 'nvt':
        return simulationEngine.Nvt(**kwargs)
    elif procedure == 'nvt-slab':
        return simulationEngine.NvtSlab(**kwargs)
    else:
        raise Exception('Unknown simulation procedure')


class PbsJob(db.Model):
    __tablename__ = 'pbs_job'
    id = NotNullColumn(Integer, primary_key=True)
    name = NotNullColumn(String(200))
    sh_file = NotNullColumn(Text)
    sh_content = Column(Text, nullable=True)
    submitted = NotNullColumn(Boolean, default=False)

    def __repr__(self):
        return '<PbsJob: %i: %s %s>' % (self.id, self.name, self.submitted)

    def submit(self):
        if jobmanager.submit(self.sh_file):
            self.submitted = True
        else:
            self.submitted = False
            log.warning('Submit PBS job failed %s' % self.sh_file)
        db.session.commit()


class Compute(db.Model):
    __tablename__ = 'compute'
    id = NotNullColumn(Integer, primary_key=True)
    web_id = NotNullColumn(Integer)
    web_user_id = NotNullColumn(Integer)
    web_ip = NotNullColumn(String(200))
    time = NotNullColumn(DateTime, default=datetime.now)
    json = NotNullColumn(Text)
    remark = NotNullColumn(Text)

    tasks = db.relationship('Task', lazy='dynamic')

    def __repr__(self):
        return '<Compute: %i: %s>' % (self.id, self.remark)

    def create_tasks(self):
        log.info('Create tasks from %s' % self)
        try:
            detail = json.loads(self.json)
            procedures = detail['procedures']
            combinations = detail['combinations']
            t_range_default = detail.get('t')
            p_range_default = detail.get('p')

            # check for procedure and prior
            for p in detail['procedures']:
                if p not in Procedure.choices:
                    raise Exception('Invalid simulation procedure: %s' % p)

                prior = Procedure.prior.get(p)
                if prior is not None and prior not in procedures:
                    procedures = [prior] + procedures  # Put prerequisite at first

            # TODO ignore duplicated smiles combinations
            # all_smiles_list = []
            for combination in combinations:
                # check for smiles and name
                smiles_list = combination['smiles']
                name_list = combination.get('names')

                # if smiles_list in all_smiles_list:
                #     raise Exception('Duplicated smiles combinations')
                # all_smiles_list.append(smiles_list)

                if len(smiles_list) > len(set(smiles_list)):
                    raise Exception('Duplicated smiles %s' % smiles_list)

                if name_list is not None:
                    if len(smiles_list) != len(name_list):
                        raise Exception('Length of names not equal to smiles')  # Provide name for no or all molecules
                    elif not re.match('^[A-Za-z0-9-]+$', ''.join(name_list)):
                        raise Exception('Only letters, numbers and - are valid for names')

                t_range = combination.get('t') or t_range_default
                p_range = combination.get('p') or p_range_default

                if t_range is not None:
                    t_min, t_max = t_range
                    if t_min > t_max:
                        raise Exception('Invalid temperature')
                else:
                    t_min, t_max = None, None

                if p_range is not None:
                    p_min, p_max = p_range
                    if p_min > p_max:
                        raise Exception('Invalid pressure')
                else:
                    p_min, p_max = None, None

                for procedure in procedures:
                    if Task.query.filter(
                            and_(Task.compute_id == self.id,
                                 Task.smiles_list == json.dumps(smiles_list),
                                 Task.procedure == procedure
                                 )
                    ).first() is not None:
                        continue

                    task = Task()
                    task.compute_id = self.id
                    task.n_components = len(smiles_list)
                    task.smiles_list = json.dumps(smiles_list)
                    task.procedure = procedure
                    if procedure in Procedure.T_RELEVANT:
                        if t_min == None or t_max == None:
                            raise Exception('Invalid temperature')
                        task.t_min = t_min
                        task.t_max = t_max
                    else:
                        task.t_min = None
                        task.t_max = None
                    if procedure in Procedure.P_RELEVANT:
                        if p_min == None or p_max == None:
                            raise Exception('Invalid pressure')
                        task.p_min = p_min
                        task.p_max = p_max
                    else:
                        task.p_min = None
                        task.p_max = None
                    if name_list is not None:
                        task.name = '_'.join(name_list) + '_' + random_string(4)
                    db.session.add(task)
            db.session.commit()
        except Exception as e:
            db.session.rollback()
            log.error('Create tasks failed %s %s' % (self, repr(e)))
            raise Exception('Create tasks failed: ' + repr(e))

    class Stage:
        SUBMITTED = 0
        BUILDING = 1
        RUNNING = 2

        text = {
            SUBMITTED: 'Submitted',
            BUILDING: 'Building...',
            RUNNING: 'Running...',
        }

    class Status:
        STARTED = 1
        DONE = 9
        FAILED = -1
        ANALYZED = 10

        text = {
            STARTED: 'Started',
            DONE: 'Done',
            FAILED: 'Failed',
            ANALYZED: 'Analyzed',
        }


class Task(db.Model):
    __tablename__ = 'task'
    id = NotNullColumn(Integer, primary_key=True)
    compute_id = NotNullColumn(Integer, ForeignKey(Compute.id))
    n_components = NotNullColumn(Integer)
    smiles_list = NotNullColumn(Text)
    n_mol_list = Column(Text, nullable=True)
    procedure = NotNullColumn(String(200))
    t_min = Column(Integer, nullable=True)
    t_max = Column(Integer, nullable=True)
    t_interval = Column(Integer, nullable=True)
    p_min = Column(Integer, nullable=True)
    p_max = Column(Integer, nullable=True)
    name = NotNullColumn(String(200), default=random_string)
    cycle = NotNullColumn(Integer, default=0)
    stage = NotNullColumn(Integer, default=Compute.Stage.SUBMITTED)
    status = NotNullColumn(Integer, default=Compute.Status.DONE)
    commands = Column(Text, nullable=True)
    remark = Column(Text, nullable=True)

    compute = db.relationship(Compute)
    jobs = db.relationship('Job', lazy='dynamic')

    def __repr__(self):
        return '<Task: %i: %s %s %s>' % (self.id, self.procedure, self.name, self.smiles_list)

    @property
    def dir(self) -> str:
        return os.path.join(Config.WORK_DIR, self.name)

    @property
    def prior_task(self):
        prior_procedure = Procedure.prior.get(self.procedure)
        if prior_procedure is None:
            return None

        # TODO Herein we suppose there is no duplicated task
        return Task.query.filter(
            and_(Task.smiles_list == self.smiles_list,
                 Task.procedure == prior_procedure,
                 Task.t_min == self.t_min,
                 Task.t_max == self.t_max,
                 Task.t_interval == self.t_interval,
                 Task.p_min == self.p_max,
                 Task.p_max == self.p_max
                 )
        ).first()

    @property
    def n_pbs_jobs(self) -> int:
        n = self.jobs.count()
        if Config.GMX_MULTIDIR:
            n = math.ceil(n / Config.GMX_MULTIDIR_NJOB)
        return n

    @property
    def n_mol_total(self) -> int:
        n_mol = 0
        for i in json.loads(self.n_mol_list):
            n_mol += i
        return n_mol

    def build(self):
        log.info('Build task %s' % self)
        try:
            cd_or_create_and_cd(self.dir)
            cd_or_create_and_cd('build')

            self.stage = Compute.Stage.BUILDING
            self.status = Compute.Status.STARTED
            db.session.commit()

            simulation = init_simulation(self.procedure)
            try:
                simulation.set_system(json.loads(self.smiles_list), n_atoms=3000)
                simulation.build()
                self.n_mol_list = json.dumps(simulation.n_mol_list)
            except:
                self.status = Compute.Status.FAILED
                db.session.commit()
                raise
            else:
                try:
                    self.prepare_jobs()
                except:
                    self.status = Compute.Status.FAILED
                    db.session.commit()
                    raise
                else:
                    self.status = Compute.Status.DONE
                    db.session.commit()
        except Exception as e:
            log.error('Build task failed %s: %s' % (self, repr(e)))

    def prepare_jobs(self):
        if not os.path.exists(os.path.join(self.dir, 'build')):
            raise Exception('Should build simulation box first')

        if self.t_min is None or self.t_max is None:
            T_list = [None]
        else:
            T_list = get_T_list_from_range(self.t_min, self.t_max, interval=self.t_interval)
        if self.p_min is None or self.p_max is None:
            P_list = [None]
        else:
            P_list = get_P_list_from_range(self.p_min, self.p_max, multiple=(2, 5))
            # remove P in the range of (1,10] bar
            P_list = list(filter(lambda x: x == 1 or x > 10, P_list))

        for p in P_list:
            for t in T_list:
                if Job.query.filter(
                        and_(Job.task_id == self.id,
                             Job.t == t,
                             Job.p == p)
                ).first() is not None:
                    continue

                job = Job()
                job.task_id = self.id
                job.t = t
                job.p = p
                job.name = '%s-%i-%i' % (self.name, t or 0, p or 0)
                db.session.add(job)

        try:
            db.session.commit()
        except:
            db.session.rollback()
            raise

        for job in self.jobs:
            job.prepare()

    def run(self, ignore_pbs_limit=False, sleep=0.2):
        log.info('Run task %s' % self)
        try:
            if self.stage != Compute.Stage.BUILDING:
                raise Exception('Incorrect stage: %s' % Compute.Stage.text[self.stage])
            elif self.status != Compute.Status.DONE:
                raise Exception('Incorrect status: %s' % Compute.Status.text[self.status])

            if not ignore_pbs_limit and jobmanager.n_running_jobs + self.n_pbs_jobs >= Config.PBS_NJOB_LIMIT:
                raise Exception('PBS_NJOB_LIMIT reached')

            os.chdir(self.dir)

            self.stage = Compute.Stage.RUNNING
            self.status = Compute.Status.STARTED
            db.session.commit()

            if not Config.GMX_MULTIDIR:
                for job in self.jobs:
                    job.run()
                    time.sleep(sleep)
            else:
                multi_dirs = [job.dir for job in self.jobs]
                multi_cmds = json.loads(self.commands)

                gmx = GMX(gmx_bin=Config.GMX_BIN)
                commands_list = gmx.generate_gpu_multidir_cmds(multi_dirs, multi_cmds,
                                                               n_parallel=Config.GMX_MULTIDIR_NJOB,
                                                               n_gpu=0,
                                                               n_thread=Config.GMX_MULTIDIR_NTHREAD)
                for i, commands in enumerate(commands_list):
                    # instead of run directly, we add a record to pbs_job
                    sh = os.path.join(self.dir, '_job.run-%i.sh' % i)
                    pbs_name = '%s-run-%i' % (self.name, i)
                    jobmanager.generate_sh(self.dir, commands, name=pbs_name, sh=sh,
                                           n_thread=Config.GMX_MULTIDIR_NTHREAD, exclusive=True)
                    # jobmanager.submit(sh)

                    pbs_job = PbsJob()
                    pbs_job.name = pbs_name
                    pbs_job.sh_file = sh
                    db.session.add(pbs_job)
                    db.session.flush()

                    # save pbs_job_id for jobs
                    for job in self.jobs[i * Config.GMX_MULTIDIR_NJOB:(i + 1) * Config.GMX_MULTIDIR_NJOB]:
                        job.pbs_job_id = pbs_job.id
                        job.status = Compute.Status.STARTED
                    db.session.commit()

                    # submit job, record if success or failed
                    pbs_job.submit()
                    time.sleep(sleep)

        except Exception as e:
            log.error('Run task failed %s %s' % (self, repr(e)))

    @property
    def ready_to_extend(self) -> bool:
        if not (self.stage == Compute.Stage.RUNNING and self.status == Compute.Status.STARTED):
            return False

        # extend the task only if all the jobs are finished (done or failed)
        for job in self.jobs:
            if job.status == Compute.Status.STARTED:
                return False

        return True

    def extend(self, ignore_pbs_limit=False, sleep=0.2):
        log.info('Extend task %s' % self)

        if not self.ready_to_extend:
            log.warning('Not ready to extend %s' % self)
            return

        extend_jobs = []
        for job in self.jobs:
            if job.need_extend:
                extend_jobs.append(job)
        if len(extend_jobs) == 0:
            log.warning('No job need to extend %s' % self)
            return

        try:
            if not Config.GMX_MULTIDIR:
                n_extend_jobs = len(extend_jobs)
            else:
                n_extend_jobs = math.ceil(len(extend_jobs) / Config.GMX_MULTIDIR_NJOB_EXTEND)

            if not ignore_pbs_limit and jobmanager.n_running_jobs + n_extend_jobs >= Config.PBS_NJOB_LIMIT:
                raise Exception('PBS_NJOB_LIMIT reached')

            self.cycle += 1

            if not Config.GMX_MULTIDIR:
                for job in extend_jobs:
                    job.extend()
                    time.sleep(sleep)
            else:
                multi_dirs = []
                multi_cmds = []
                simulation = init_simulation(self.procedure)
                for job in extend_jobs:
                    os.chdir(job.dir)
                    multi_dirs.append(job.dir)
                    multi_cmds = simulation.extend(jobname='%s-%i' % (job.name, job.cycle + 1),
                                                   sh='_job.extend-%i.sh' % (job.cycle + 1))

                commands_list = simulation.gmx.generate_gpu_multidir_cmds(multi_dirs, multi_cmds,
                                                                          n_parallel=Config.GMX_MULTIDIR_NJOB_EXTEND,
                                                                          n_gpu=0,
                                                                          n_thread=Config.GMX_MULTIDIR_NTHREAD)

                os.chdir(self.dir)
                for i, commands in enumerate(commands_list):
                    # instead of run directly, we add a record to pbs_job
                    sh = os.path.join(self.dir, '_job.extend-%i-%i.sh' % (self.cycle, i))
                    pbs_name = '%s-extend-%i-%i' % (self.name, self.cycle, i)
                    jobmanager.generate_sh(self.dir, commands, name=pbs_name, sh=sh,
                                           n_thread=Config.GMX_MULTIDIR_NTHREAD, exclusive=True)
                    # jobmanager.submit(sh)

                    pbs_job = PbsJob()
                    pbs_job.name = pbs_name
                    pbs_job.sh_file = sh
                    db.session.add(pbs_job)
                    db.session.flush()

                    # save pbs_job_id for jobs
                    for job in extend_jobs[i * Config.GMX_MULTIDIR_NJOB_EXTEND:(i + 1) * Config.GMX_MULTIDIR_NJOB_EXTEND]:
                        job.pbs_job_id = pbs_job.id
                        job.cycle += 1
                        job.status = Compute.Status.STARTED
                    db.session.commit()

                    # submit job, record if success or failed
                    pbs_job.submit()
                    time.sleep(sleep)

        except Exception as e:
            log.error('Extend task failed %s %s' % (self, repr(e)))

    def check_finished(self):
        """
        check if all jobs in this tasks are finished
        if finished, analyze the job
        """
        log.info('Check task status %s' % self)
        try:
            for job in self.jobs:
                try:
                    job.check_finished()
                except Exception as e:
                    log.error('Check job status failed %s %s' % (job, repr(e)))

                if job.status == Compute.Status.DONE:
                    try:
                        log.info('Analyze job %s' % job)
                        job.analyze()
                    except Exception as e:
                        log.error('Analyze job failed %s %s' % (job, repr(e)))

            # Set status as DONE only if all jobs are converged
            for job in self.jobs:
                if not job.converged:
                    break
            else:
                self.status = Compute.Status.DONE
                db.session.commit()
        except Exception as e:
            log.error('Check task status failed %s %s' % (self, repr(e)))

    def remove(self):
        for job in self.jobs:
            job.remove()
        try:
            shutil.rmtree(self.dir)
        except:
            log.warning('Remove task %s Cannot remove folder: %s' % (self, self.dir))

        db.session.delete(self)
        db.session.commit()

    def get_isothermal_result(self, T=298) -> ([int], [float], [float]):
        P_list = []
        result_list = {}
        stderr_list = {}
        jobs = self.jobs.filter(Job.t == T)
        for job in jobs:
            if job.converged:
                P_list.append(job.p)
                results = json.loads(job.result)
                for k, v in results.items():
                    if not k in result_list:
                        result_list[k] = []
                        stderr_list[k] = []
                    if type(v) == list:
                        result_list[k].append(v[0])
                        stderr_list[k].append(v[1])
                    else:
                        result_list[k].append(v)

                    if k == 'e_inter':
                        if not 'hvap' in result_list:
                            result_list['hvap'] = []
                            stderr_list['hvap'] = []
                        result_list['hvap'].append(8.314 * job.t / 1000 - v[0] / self.n_mol_total)
                        stderr_list['hvap'].append(v[1] / self.n_mol_total)

        return P_list, result_list, stderr_list

    def get_isobaric_result(self, P=1) -> ([int], [float], [float]):
        T_list = []
        result_list = {}
        stderr_list = {}
        jobs = self.jobs.filter(Job.p == P)
        for job in jobs:
            if job.converged:
                T_list.append(job.t)
                results = json.loads(job.result)
                for k, v in results.items():
                    if not k in result_list:
                        result_list[k] = []
                        stderr_list[k] = []
                    if type(v) == list:
                        result_list[k].append(v[0])
                        stderr_list[k].append(v[1])
                    else:
                        result_list[k].append(v)

                    if k == 'e_inter':
                        if not 'hvap' in result_list:
                            result_list['hvap'] = []
                            stderr_list['hvap'] = []
                        result_list['hvap'].append(8.314 * job.t / 1000 - v[0] / self.n_mol_total)
                        stderr_list['hvap'].append(v[1] / self.n_mol_total)

        return T_list, result_list, stderr_list


class Job(db.Model):
    __tablename__ = 'job'
    id = NotNullColumn(Integer, primary_key=True)
    task_id = NotNullColumn(Integer, ForeignKey(Task.id))
    t = Column(Integer, nullable=True)
    p = Column(Integer, nullable=True)
    time = NotNullColumn(DateTime, default=datetime.now)
    name = NotNullColumn(String(200), default=random_string)
    cycle = NotNullColumn(Integer, default=0)
    status = NotNullColumn(Integer, default=Compute.Status.STARTED)
    converged = NotNullColumn(Boolean, default=False)
    result = Column(Text, nullable=True)
    pbs_job_id = Column(Integer, ForeignKey(PbsJob.id), nullable=True)

    task = db.relationship(Task)
    pbs_job = db.relationship(PbsJob)

    def __repr__(self):
        return '<Job: %i: %s %i>' % (self.id, self.name, self.cycle)

    @property
    def dir(self) -> str:
        dir_name = '%i-%i' % (self.t or 0, self.p or 0)
        return os.path.join(self.task.dir, dir_name)

    @property
    def prior_job(self):
        prior_task = self.task.prior_task
        if prior_task is None:
            return None

        for job in prior_task.jobs:
            if job.t == self.t and job.p == self.p:
                return job

    @property
    def pbs_job_done(self) -> bool:
        if self.pbs_job_id is None:
            return False
        elif not self.pbs_job.submitted:
            return False
        else:
            return not jobmanager.is_running(self.pbs_job.name)

    @property
    def need_extend(self) -> bool:
        return self.status == Compute.Status.ANALYZED and not self.converged

    def prepare(self):
        prior_job = self.prior_job
        if prior_job is None:
            prior_job_dir = None
        else:
            if not prior_job.converged:
                log.warning('Prepare job failed %s Prior job not converged' % repr(self))
                return
            else:
                prior_job_dir = prior_job.dir

        cd_or_create_and_cd(self.dir)
        simulation = init_simulation(self.task.procedure)
        commands = simulation.prepare(model_dir='../build', T=self.t, P=self.p, jobname=self.name,
                                      prior_job_dir=prior_job_dir, drde=True)  # Temperature dependent parameters

        # all jobs under one task share same commands, save this for GMX -multidir simulation
        if self.task.commands is None:
            self.task.commands = json.dumps(commands)
            db.session.commit()

    def run(self):
        try:
            os.chdir(self.dir)
        except:
            raise Exception('Should prepare job first')

        # instead of run directly, we add a record to pbs_job

        # simulation = init_simulation(self.task.procedure)
        # simulation.run()

        pbs_job = PbsJob()
        pbs_job.name = self.name
        pbs_job.sh_file = os.path.join(self.dir, jobmanager.sh)
        db.session.add(pbs_job)
        db.session.flush()

        self.pbs_job_id = pbs_job.id
        self.status = Compute.Status.STARTED
        db.session.commit()

        # submit job, record if success or failed
        pbs_job.submit()

    def extend(self):
        if not self.need_extend:
            log.warning('Will not extend %s Status: %s' % (self, Compute.Status.text[self.status]))
            return

        os.chdir(self.dir)

        pbs_name = '%s-%i' % (self.name, self.cycle + 1)
        sh = os.path.join(self.dir, '_job.extend-%i.sh' % (self.cycle + 1))

        # instead of run directly, we add a record to pbs_job

        simulation = init_simulation(self.task.procedure)
        simulation.extend(jobname=pbs_name, sh=sh)
        # simulation.run()

        pbs_job = PbsJob()
        pbs_job.name = pbs_name
        pbs_job.sh_file = sh
        db.session.add(pbs_job)
        db.session.flush()

        self.pbs_job_id = pbs_job.id
        self.cycle += 1
        self.status = Compute.Status.STARTED
        db.session.commit()

        # submit job, record if success or failed
        pbs_job.submit()

    def check_finished(self) -> bool:
        if self.status in (Compute.Status.DONE, Compute.Status.FAILED, Compute.Status.ANALYZED):
            return True
        if not self.pbs_job_done:
            return False

        os.chdir(self.dir)

        simulation = init_simulation(self.task.procedure)
        if simulation.check_finished():
            self.status = Compute.Status.DONE
        else:
            self.status = Compute.Status.FAILED
            log.error('Job failed %s' % self)
        db.session.commit()
        return True

    def analyze(self, rerun=False):
        if self.status == Compute.Status.ANALYZED and not rerun:
            log.warning('Will not analyze %s Already analyzed' % self)
            return
        if self.status == Compute.Status.STARTED:
            log.warning('Will not analyze %s Job is still running' % self)
            return
        if self.status == Compute.Status.FAILED:
            log.warning('Will not analyze %s Job failed' % self)
            return

        os.chdir(self.dir)

        simulation = init_simulation(self.task.procedure)
        dirs = [self.dir]
        try:
            result = simulation.analyze(dirs)
        except:
            self.status = Compute.Status.FAILED
            db.session.commit()
            raise
        else:
            self.status = Compute.Status.ANALYZED
            if result is None:
                self.converged = False
                self.result = None
            else:
                self.converged = True
                self.result = json.dumps(result)
            db.session.commit()

    def remove(self):
        try:
            jobmanager.kill_job(self.pbs_name)
        except Exception as e:
            log.warning('Remove job %s Cannot kill PBS job: %s' % (self, repr(e)))

        try:
            shutil.rmtree(self.dir)
        except:
            log.warning('Remove job %s Cannot remove folder: %s' % (self, self.dir))

        db.session.delete(self)
        db.session.commit()
