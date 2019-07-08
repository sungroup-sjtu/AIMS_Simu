# Molecule Simulation Database -- Server
This project is performs high-throughput force field simulation and data processing  
This project depends on `AIMS_Tools`  
Two main scripts accomplish the jobs **submit.py** and **monitor.py**.

## To setup on a Linux server

* clone AIMS_Simu and AIMS_Tools

* Modify `config.py`  
  a) Set paths of MS_TOOLS, WORK_DIR, PACKMOL, DFF and DFF database.  
  b) Select job queue (fast, gtx, cpu on L86)  
  c) Set GROMACS executable accordingly.
  
  ```
  MS_TOOLS_DIR = os.path.join(CWD, '../AIMS_Tools')
  WORK_DIR = '/share/md1400/_MSDServer/'
  PACKMOL_BIN = '/share/apps/tools/packmol'
  DFF_ROOT = '/home/gongzheng/apps/DFF/Developing'
  DFF_TABLE = 'MGI'
  ```

* Run to prepare a high-throughput computation:  
    `run/submit.py [npt,nvt-slab] mols/example.txt 'comments'`
* Run to perform automatic job flow 
    `run/monitor.py [npt, nvt-slab]` 

  ```
  ./submit.py npt mols/example.txt 'test computation'
  ./monitor.py npt 
  ```

See our publication for details
`Predicting Thermodynamic Properties of Alkanes by High-throughput Force Field Simulation and Machine Learning`
https://doi.org/10.1021/acs.jcim.8b00407
