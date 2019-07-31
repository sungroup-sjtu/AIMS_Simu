# Molecule Simulation Database -- Server
This project is performs high-throughput force field simulation and data processing  
This project depends on `AIMS_Tools`  
Two main scripts accomplish the jobs **submit.py** and **monitor.py**.

## To setup on a Linux server

1. Clone AIMS_Simu and AIMS_Tools:  
    `$git clone https://github.com/sungroup-sjtu/AIMS_Simu'

2. Modify `config.py`  
  a) Set paths of MS_TOOLS, WORK_DIR, PACKMOL, DFF and DFF database.  
  b) Select the job queue (fast, gtx, cpu on Cluster 86)  
  c) Set GROMACS executable  

    For example:
    ```
    MS_TOOLS_DIR = os.path.join(CWD, '../AIMS_Tools')
    WORK_DIR = '/share/md1400/_MSDServer/'
    PACKMOL_BIN = '/share/apps/tools/packmol'
    DFF_ROOT = '/home/gongzheng/apps/DFF/Developing'
    DFF_TABLE = 'MGI'
    ...
    PBS_ARGS = ('gtx', 32, 2, 16)  # partition, cpu, gpu, cpu_request
    ...
    GMX_BIN = '/share/apps/gromacs/2016.6/bin/gmx_serial'
    GMX_MDRUN = 'gmx_gpu mdrun'
    # GMX_MDRUN= 'gmx_fast mdrun'
    GMX_MULTI = True
    GMX_MULTI_NJOB = 8  # Use -multidir function of GROMACS. For Npt simulation, set it to 8. For NvtSlab simulation, 4 is better
    GMX_MULTI_NOMP = None  # Set the OpenMP threads. When set to None, use only one node and the best number of threads is automatically determined

    ```

3. To set up a high-throughput multi-task computation. Prepare a list 
of molecules in /mols/<fname> and then:  
    `run/submit.py [npt,nvt-slab] mols/example.txt 'comments'`
    
4. To run the calculations:   
    `run/monitor.py [npt, nvt-slab]`

    For example:
    ```
    cd run
    ./submit.py npt mols/example.txt 'testing'
    ./monitor.py npt 
    ```
5. The results are saved in the WORK_DIR


For more information, see our publication: "Predicting 
Thermodynamic Properties of Alkanes by High-throughput 
Force Field Simulation and Machine Learning", 
https://doi.org/10.1021/acs.jcim.8b00407
