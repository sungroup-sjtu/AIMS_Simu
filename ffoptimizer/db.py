from sqlalchemy.ext.declarative import declarative_base
from sqlalchemy.orm import sessionmaker
from sqlalchemy import create_engine

class DB:
    DBFILE = 'ffoptimizer.db'
    Base = declarative_base()
    engine = None
    session = None

    @staticmethod
    def conn():
        url = 'sqlite:///' + DB.DBFILE
        try:
            DB.engine = create_engine(url, echo=False)
            Session = sessionmaker(DB.engine)
            DB.session = Session()
        except Exception as e:
            print(str(e))
            return False

        try:
            DB.Base.metadata.create_all(DB.engine)
        except Exception as e:
            print(str(e))
            return False
        else:
            return True

    @staticmethod
    def close():
        DB.session.close()
        DB.engine.dispose()
