# CFD-based WFLO

Optimization framework for designing wind farms using computational fluid dynamics.


A description of the framework can be found in the following papers:

Enrico G.A. Antonini, David A. Romero, Cristina H. Amon. "Continuous adjoint formulation for wind farm layout optimization: A 2D implementation", Applied Energy, 2018;

Enrico G.A. Antonini, David A. Romero, Cristina H. Amon. "Optimal design of wind farms in complex terrains using computational fluid dynamics and adjoint methods", Applied Energy, 2020.


OpenFOAM 6 (https://openfoam.org/release/6/), NLopt (https://nlopt.readthedocs.io/en/latest/), Python and Open MPI need to be installed as they are fundamental components of the framework.



TO COMPILE THE CODE

cd ./Codes;
make;

It will generate the following files:
- CDGradientCalculation (to calulate the gradient using a central difference scheme - not needed for the optimization using adjoint methods)
- GradientCalculation (to calulate the gradient using the adjoint method)
- PowerCalculation (to calulate the annual energy production)
- SLSQP (optimizer)

(If you need to manually specify OpenFOAM or MPI binary file paths, go to line  

TO RUN THE CODE

mkdir TestCase;
cp ./Codes/GradientCalculation ./Codes/PowerCalculation ./Codes/SLSQP ./Codes/RemapTerrain.py ./Codes/RemapTurbines.py ./TestCase/.;
cd ./TestCase;
./SLSQP;
