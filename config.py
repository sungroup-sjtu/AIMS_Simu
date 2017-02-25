import os


class Config:
    DB_NAME = os.path.dirname(os.path.abspath(__file__)) + os.sep + 'msdserver.sqlite'

    SQLALCHEMY_DATABASE_URI = 'sqlite:///%s' % DB_NAME
    SQLALCHEMY_COMMIT_ON_TEARDOWN = True
    SQLALCHEMY_TRACK_MODIFICATIONS = True

    WORK_DIR = '/tmp/MSDataServer/'
    USER = 'zheng'
