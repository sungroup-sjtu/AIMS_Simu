class Node():
    def __init__(self, name, queue, nprocs, free_procs, state):
        self.name = name
        self.queue = queue
        self.nprocs = nprocs
        self.free_procs = free_procs
        self.state = state
