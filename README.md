# Molecule Simulation Database -- Server
This project performs high-throughput force field simulation and data processing  
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
1. Fitting the simulation data at different temperature and pressure. So that properties and derivatives at arbitrary T or P can be obtained. Normally this should be prior to any other analyzing  
  `./scripts-post/post-process.py [npt, nvt-slab]`
2. Dump the simulation data from sqlite database to `csv` file, which can be uploaded into `AIMS_Web` database  
  `./scripts/dummp-data-npt.py`  
  `./scripts/dummp-data-nvt-slab.py`
2. Compare with NIST Experimental data  
   * Make sure that `nist.sqlite` exists in `database` folder  
   * Prepare a file which lists the SMILES of molecules you want to compare. An example is given as `smiles-nist.txt`
   * Run following script to compare simulation and expt data and plot the results  
     `./scripts-post/compare-nist.py [npt, nvt-slab] smiles-nist.txt`
  