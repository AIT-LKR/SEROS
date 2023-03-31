[![DOI](https://zenodo.org/badge/620307347.svg)](https://zenodo.org/badge/latestdoi/620307347)
# SEROS
This tool runs a Cellular Automaton (CA) combining the Lattice-Boltzmann Method (LBM) for fluid flow simulation with an heuristic reshaping algorithm mimicking SEdimentation and eROSion (hence the name 'SEROS').
The following animations show the shear stress (left), the density (center) and the velocity (right) fields during the reshaping of a 2D elbow channel (taken from case '2D_L_Re40_x10_s15_r0.02_D'):

<img src="./examples/shearStress.gif" alt="shear stress field" width="250"/> <img src="./examples/density.gif" alt="topology" width="250"/> <img src="./examples/velocity.gif" alt="velocity field" width="250"/>

## Acknowledgement
- SEROS is utilising a [modified version of the Palabos library](https://github.com/AIT-LKR/palabos), a framework for general-purpose computational fluid dynamics (CFD), with a kernel based on the lattice Boltzmann (LB) method. More information can be found here https://palabos.unige.ch/ and the most recent version of the code is hosted by the developers here: https://gitlab.com/unigespc/palabos.
- The reshaping constraint for a constant diameter uses [IsoSurfaceExtraction](https://github.com/AIT-LKR/IsoSurfaceExtraction) a fork of https://github.com/mkazhdan/IsoSurfaceExtraction for calculation of the surface area of the channel.
- For an experimental feature to load 2D domains from a ppm, the code https://github.com/fmenozzi/easyppm was also added to the "externalLibraries". This is no crucial part to use SEROS but it is necessary to be in the source code for compilation.

## Funding
The work presented here was conducted in the context of the Clean Sky 2 project COMBO3D. This project has received funding from the Clean Sky 2 Joint Undertaking under the European Union’s Horizon 2020 research and innovation programme under grant agreement No [831851](https://doi.org/10.3030/831851). This publication reflects only the author’s views and the European Union is not liable for any use that may be made of the information contained therein.

## Setup
### Dependencies
install mpicxx:
```cmd
sudo apt install make python2 g++
sudo apt install openmpi-bin libopenmpi-dev
```
install Gnu Scientific Library:
```cmd
sudo apt install libgsl-dev
```

### Compilation
If you used a compressed version of this software, all necessary files might be present.
If you use a cloned repo, you might need to use:
```cmd
git submodule update --init --recursive --remote
```
to load the files located in submodules.
You also need to compile an additional executable for the reshaping constraint to work properly which is located in:
```cmd
cd /palabos/externalLibraries/IsoSurface
make
```
The actual solver-reshaping executable is compiled with:
```cmd
cd pipe3D
make
```

## Usage
### Execution
Now, you can directly start a simulation and reshaping of a 2-dimensional elbow shaped channel case with:
```cmd
./seros3D config2D_L.xml
```
The Reynolds number (Re) and resolution (x) can also be specified as arguments, for example Re=40 and x=10 nodes across inlet:
```cmd
./seros3D config2D_L.xml 40 10
```

#### Multi-processing
For 3-dimensional cases with high resolution, it is recommended to use multiple processors. To start the 3D elbow case with 4 processes, run:
```cmd
mpirun -np 4 ./seros3D config3D_L.xml 40 40
```

### Output
- Logs
  - The executable will print some numerical and physical details about the case to the terminal and save it in the file 'plb.log' and 'plb.error'.
  - In the file 'bestSeros', the algorithm will store the reshaping iteration number at which the pressure drop was the smallest.
  - 'extractorLog' contains the output of the 'IsoSurface' utility used to calculate the surface area of the channel, which is used for maintaining the initial diameter during reshaping.
    The surfaces are stored in the files 'IsoSurface_init.ply' (initial geometry), 'IsoSurface.ply' (latest reshaping iteration) 'IsoSurface.ply' (iteration yielding smallest pressure drop during current execution)
- Data
  - Visualisation
    - For quick assesment of the fluid flow during the reshaping, images are written to the disk in the PPM format.
    - For in-depth post processing of the flow field, vtk files of the flow  ('vtk#####.vti) and the topology ('vtkInt#####.vti) are also stored:
      - 'vtk00000000_init.vti' for the time step at which a sufficiently steady flow was reached,
      - 'vtk000000-1_seros-1.vti' for the topology which yielded the smallest pressure drop,
      - and at specified iteration intervalls (e.g. 'vtk00303000_seros100.vti' after 100 reshaping iterations).
  - Numerical data:
    - The geometry (topology) is saved in files starting with the name 'topo' both in binary and ASCII format. They can be used for debugging and resuming the simulation. 
    - The flow field for an associated topology is saved in the file starting with 'popu' containing the Lattice-Boltzmann populations for each cell. For each simulation, the flow field has to be initialised 
  - Fluid flow data:
    - measurementsFlow_initialize contains fluid flow data recorded at particular time steps during the start of the flow simulation until it has reached a steady state.  measurementsFlow0 Then records the fluid flow data during the reshaping.
    - measurementsSeros0 contains statistical fluid flow data and geometrical details for each reshaping iteration. 
    - If the execution is run again for the same case, it will create new files to contain the new data (e.g. measurementsFlow1 and measurementsSeros1)
    - 'analysis.gnu' can be used to visualise the flow data during and after the execution using gnuplot:
      ```cmd
      gnuplot analysis.gnu
      ```
    - The distributions of shear stress values in the wall near region are also stored (in the directories 'distribution'), and can be converted to a histogram (using 'postprocessing/histogram.py'). However, no further use has been made of these files yet.

### Post processing
#### Generate animated images from PPMs
Copy and execute this script within a case directory:
```cmd
bash generateAnimations.sh
```
It will try to generate .gif and .avi files and put them in a new folder called 'animations' within the case directory.
This process might require additional packages to be installed.

#### Comparing multiple cases
A python script to evaluate all cases contained in the 'pipe3D' directory in batch can be executed with
```cmd
./pipe3D/postprocessing/python/analysis.py ./pipe3D/
```
This will create a file 'rhobarInletAvg_analysis' like this example:
```cmd
simulation best_seros init_inlet end_inlet space init_width end_width
2D_L_Re40_x10_s15_r0.02_D 71.0 0.00985531 0.00904744 0 1.0 1.0052
```
The values are derived by the python script from the file 'measurementsSeros' in the case directories and have the following meaning:
- simulation (name of case)
- best_seros (iteration with min pressure drop)
- init_inlet (initial pressure at inlet)
- end_inlet (pressure at inlet for best_seros)
- space (placeholder with value 0, can be ignored)
- init_width (dimensionless of inlet)
- end_width (for best_seros)

