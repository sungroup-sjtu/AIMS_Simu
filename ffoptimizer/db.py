from sqlalchemy.orm import sessionmaker
from sqlalchemy import create_engine

from ffoptimizer.target import metadata


class DB():
    def __init__(self, dbfile):
        self.dbfile = dbfile
        self.engine = None
        self.session = None

    def conn(self):
        url = 'sqlite:///' + self.dbfile
        try:
            self.engine = create_engine(url, echo=False)
            Session = sessionmaker(self.engine)
            self.session = Session()
        except Exception as e:
            print(str(e))
            return False

        try:
            metadata.create_all(self.engine)
        except Exception as e:
            print(str(e))
            return False
        else:
            return True

    def close(self):
        self.session.close()
        self.engine.dispose()
