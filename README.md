# MRST Simulations

My various reservoir simulations done in Matlab Reservoir Simulation Toolbox (MRST). Can be run in Octave.

Simple 3D homogeneous/heterogeneous models:

1. Single-phase incompressible fluid (water) simulation (ALL B.C.)
2. Single-phase compressible fluid (oil) simulation - constant viscosity over pressure (NO FLOW B.C.)
3. Single-phase compressible fluid (gas) simulation - pressure-dependent viscosity (NO FLOW B.C.)
4. Single-phase compressible fluid (polymer) simulation - Non-Newtonian fluid (with 2 numerical methods: cell-based, or face-based) (NO FLOW B.C.)
5. Single-phase compressible fluid (oil) simulation - constant viscosity over pressure (NO FLOW B.C.)
6. Single-phase compressible fluid simulation with thermal effect (NO FLOW B.C.)

> These scripts are adapted from MRST tutorial codes, but I have modified it.

**Note:**

* I added some Matlab functionalities used in the scripts that Octave currently doesn't have, such as `deval`. See in the `./modules/nuwara` folder. 
* I experienced that `plotyy` (an Octave function) doesn't work properly. `plotyy` is used e.g. in Simulation 2 to 6 (above). To fix this, this is my way round (in the command line): 
  * Add path: `addpath 'C:\Octave\Octave-5.2.0\mingw64\share\octave\5.2.0\m\plot\draw'`
  * Run: `plotyy`
