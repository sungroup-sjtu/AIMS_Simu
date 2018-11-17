# Molecule Simulation Database -- Server
This project is used for performing high-throughput force field simulation and data post-processing

* Main functions are located at `app/models.py`
* Modify `config.py` to make it works and ensure the performance
* Run `run/submit.py` to submit a high-throughput computation
* Run `run/monitor.py` to perform automatic model building, running, analyzing and extending

See our publication for details
`Predicting Thermodynamic Properties of Alkanes by High-throughput Force Field Simulation and Machine Learning`
https://doi.org/10.1021/acs.jcim.8b00407

This project relies on `ms-tools`
https://github.com/z-gong/ms-tools

The maching learning part is located at
https://github.com/z-gong/mdlearn
