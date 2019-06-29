# Molecule Simulation Database -- Server
This project performs high-throughput force field simulation and data processing  
This project depends on `AIMS_Tools`  
Two main scripts accomplish the jobs **submit.py** and **monitor.py**.

## To setup on a Linux server

1. Clone AIMS_Simu and AIMS_Tools:  
    It is suggested to put them in same folder.
    ```
    git clone https://github.com/sungroup-sjtu/AIMS_Simu  
    git clone https://github.com/sungroup-sjtu/AIMS_Tools
    ```
2. Modify `config.py`  
    a) Set paths of MS_TOOLS, WORK_DIR, PACKMOL, DFF and DFF database.  
    b) Select the job queue (fast, gtx, cpu on Cluster 86)  
    c) Set GROMACS executable  

    For example:
    ```
    force field setting(no default setting):
        DFF_TABLE = 'MGI' # 'IL'
    paths setting:
        MS_TOOLS_DIR = os.path.join(CWD, '..', 'AIMS_Tools') # AIMS_Simu and AIMS_Tools in same folder
        WORK_DIR = os.path.join(CWD, 'SimulationData') # all simulation data is saved in a new folder AIMS_Simu/SimulationData
        PACKMOL_BIN = '/share/apps/tools/packmol'
        DFF_ROOT = '/share/workspace/xiangyan/src/DFF/Developing' # simulation paramters come from this folder
    PBS settings:
        ...
        PBS_ARGS = ('gtx', 32, 2, 16)  # partition, cpu, gpu, cpu_request
        ...
        GMX_BIN = '/share/apps/gromacs/2018.6/bin/gmx_serial'
        GMX_MDRUN = 'gmx_gpu mdrun'
        # GMX_MDRUN= 'gmx_fast mdrun'
        GMX_MULTI = True
        GMX_MULTI_NJOB = 8  # Use -multidir function of GROMACS. For Npt simulation, set it to 8. For NvtSlab simulation, 4 is better
        GMX_MULTI_NOMP = None  # Set the OpenMP threads. When set to None, use only one node and the best number of threads is automatically determined
    simulation details settings (default setting is OK):
        NATOMS = 3000 # least number of atoms build in simulation box.
        NMOLS = 120 # least number of molecules build in simulation box.
        LJ96 = False # using LJ 9-6 non-bonded potential
        DIFF_GK = False # using green-kubo method to calculate the diffusion constant. (Expensive, not suggest)
        DEBUG = False # if true: do not delete the trajectory file in analyze process.
        class NvtMultiConfig(Config, SunRunConfig, SunExtendConfig, SunBugFixConfig):
            REPEAT_NUMBER = 80 # set the number of parallel simulation for nvt-multi
    ```

3. To set up a high-throughput computation. Prepare a list of molecules in /mols/<fname> and then:  
    the example.txt contains 5 columns: name SMILES molecular_ratio t_list p_list  
    more than 5 points for t_list and p_list is needed, otherwise some analysis and dumps scripts will not work.  
    `run/submit.py -p [npt,nvt-slab] -i mols/example.txt -r 'comments' -tp assigned`
    
4. To run the calculations:   
    `run/monitor.py -p [npt, nvt-slab]`
    
    Available procedures: npt, nvt-slab, ppm(npt), nvt-multi(npt).  
    ppm(npt) means npt is prerequisite of ppm.  
    Example for npt procedure:
    ```
    cd run
    ./submit.py -p npt -i mols/example.txt -r testing -tp assigned
    ./monitor.py -p npt 
    ```

5. The results are saved in the WORK_DIR. In default, WORK_DIR=AIMS_Simu/SimulationData.

## QM calculation
QM calculation for heat capacity is performed by `run-cv.py` script

1. Prepare QM files. This will check database to remove duplicated molecules from `example.txt` and process the molecule name. A file named `_cv_prepared.txt` will be generated and used for following steps  
  `./run-cv.py prepare mols/example.txt`  
2. Generate Gauss input files and submit to PBS job manager  
  `./run-cv.py cv _cv_prepared.txt`  
3. Analyze Gauss results. The results will be saved in a file named `_cv.log`  
  `./run-cv.py get-cv _cv_prepared.txt`
4. Save results into database  
  `./run-cv.py save-db`

For more information, see our publication: "Predicting 
Thermodynamic Properties of Alkanes by High-throughput 
Force Field Simulation and Machine Learning", 
https://doi.org/10.1021/acs.jcim.8b00407

## Scripts for data post-processing and analyzing
Several scripts are provided for post-processing and analysing the simulation data. They are located at `scripts` and `scripts-post`
1. Fitting the simulation data at different temperature and pressure. So that properties and derivatives at arbitrary T or P can be obtained.  
   **This should be performed prior to any other analyzing**  
  `./scripts/post-process.py -p [npt, nvt-slab] -o True`
2. Remark molecules containing specific groups (e.g. halide, cyclo-ester) as `bad` molecules, which will not be dumped in following steps  
  `./scripts/remark.py [npt, nvt-slab]`
3. Dump the molecules from sqlite database to `mols.csv` file. The category should be specified, which is necessary for uploading to `AIMS_Web` database  
  `./scripts/dummp-mols.py [small molecule, ionic liquid, ...]`  
4. Dump the simulation data from sqlite database to `csv` file, which can be uploaded into `AIMS_Web` database  
  `./scripts/dummp-data-npt.py`  
  `./scripts/dummp-data-nvt-slab.py`
5. You select specific class of molecules in following analysis based on force field atom type, by modify the app/selection.py.
6. Compare with NIST Experimental data  
   * Make sure that `nist.sqlite` exists in `database` folder  
   * Run following script to compare simulation and expt data and plot the results  
        ```
        cd scrips-post
        python3 compare_detail.py -p npt -t nist --selection False
        python3 compare-nist.py -p npt --selection False
        ```
7. Compare with ILTHERMO Experimental data 
   * Make sure that `ilthermo.sqlite` exists in `database` folder  
   * Run following script to compare simulation and expt data and plot the results  
        ```
        cd scrips-post
        python3 compare_detail.py -p npt -t ilthermo --selection False
        python3 compare-nist-il.py -p npt --selection False
        ```