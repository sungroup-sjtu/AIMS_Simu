from . import main

@main.route('/')
def index():
    return 'This is Molecule Simulation Data -- Server'
