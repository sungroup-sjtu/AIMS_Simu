# Molecule Simulation Database -- Server
This project is used for performing high-throughput force field simulation and data post-processing

This project relies on `ms-tools`

## Steps

* Specify paths of required packages in `config.py`  
  Be careful about GROMACS and queue to ensure the performance  
```
  MS_TOOLS_DIR = os.path.join(CWD, '../ms-tools')
  WORK_DIR = '/share/md1400/_MSDServer/'
  PACKMOL_BIN = '/share/apps/tools/packmol'
  DFF_ROOT = '/home/gongzheng/apps/DFF/Developing'
  DFF_TABLE = 'MGI'
```
* Run `run/submit.py` to submit a high-throughput computation  
```
  cd run
  ./submit.py npt mols/example.txt 'test computation'
```
* Run `run/monitor.py` to perform automatic model building, running, analyzing and extending

See our publication for details
`Predicting Thermodynamic Properties of Alkanes by High-throughput Force Field Simulation and Machine Learning`
https://doi.org/10.1021/acs.jcim.8b00407

This project relies on `ms-tools`
https://github.com/sungroup-sjtu/ms-tools

The maching learning part is located at
https://github.com/sungroup-sjtu/mdlearn
