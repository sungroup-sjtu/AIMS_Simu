from .gmx import GmxSimulation


class Npt(GmxSimulation):
    def __init__(self, **kwargs):
        super().__init__(**kwargs)
        self.procedure = 'nvt'
        self.requirement = []
        self.logs = ['nvt.log', 'cv.log']
