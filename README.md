# Molecule Simulation Database -- Server
This project is performs high-throughput force field simulation and data processing  
This project depends on `AIMS_Tools`  
Two main scripts accomplish the jobs **submit.py** and **monitor.py**.

## To setup on a Linux server

1. Clone AIMS_Simu and AIMS_Tools

2. Modify `config.py`  
  a) Set paths of MS_TOOLS, WORK_DIR, PACKMOL, DFF and DFF database.  
  b) Select a job queue (fast, gtx, cpu on L86)  
  c) Set GROMACS executable  

    For example:
    ```
    MS_TOOLS_DIR = os.path.join(CWD, '../AIMS_Tools')
    WORK_DIR = '/share/md1400/_MSDServer/'
    PACKMOL_BIN = '/share/apps/tools/packmol'
    DFF_ROOT = '/home/gongzheng/apps/DFF/Developing'
    DFF_TABLE = 'MGI'
    ```

3. To prepare a high-throughput multi-task computation:  
    `run/submit.py [npt,nvt-slab] mols/example.txt 'comments'`
4. To run the calculations:   
    `run/monitor.py [npt, nvt-slab]`

    For example:
    ```
    cd run
    ./submit.py npt mols/example.txt 'test computation'
    ./monitor.py npt 
    ```

For more information, see our publication: "Predicting 
Thermodynamic Properties of Alkanes by High-throughput 
Force Field Simulation and Machine Learning", 
https://doi.org/10.1021/acs.jcim.8b00407
