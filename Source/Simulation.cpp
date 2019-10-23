#include "Simulation.h"
#include <fstream>
#include <sstream>
#include <string>
#include <vector>
#include <math.h>
#include <stdlib.h>
#include <stdio.h>
#include <iostream>
#include <sys/stat.h>
#include <sys/types.h>

Simulation::Simulation() {
}

void Simulation::control(int nIter, int turb, std::string directory) {
	std::string directoryS = directory+"/system/controlDict";
	std::ofstream ofs(directoryS.c_str());
	ofs << "FoamFile{version 2.0; format ascii; class dictionary; location \"system\"; object controlDict;}\n";
	ofs << "application ";
	if (turb > 0) {
		ofs << "simpleFoamS";
	}
	else {
		ofs << "icoFoamS";
	}
	ofs << "; startForm startTime; startTime 0; stopAt endTime; endTime "<<nIter<<"; deltaT 1;\n";
	ofs << "writeControl timeStep; writeInterval "<<nIter<<";\n";
	ofs << "purgeWrite 0; writeFormat ascii; writePrecision 8; writeCompression off;\n";
	ofs << "timeFormat general; timePrecision 8; runTimeModifiable true; libs (\"libatmosphericModels.so\");\n";
	ofs.close();
}

void Simulation::a_control(int nIter, int turb, std::string directory) {
	std::string directoryS = directory+"/system/controlDict";
	std::ofstream ofs(directoryS.c_str());
	ofs << "FoamFile{version 2.0; format ascii; class dictionary; location \"system\"; object controlDict;}\n";
	ofs << "application ";
	if (turb == 1) {
		ofs << "adjointSimpleFoamT";
	}
	else if (turb == 2) {
		ofs << "adjointSimpleFoam";
	}
	else {
		ofs << "adjointIcoFoam";
	}
	ofs << "; startForm startTime; startTime 0; stopAt endTime; endTime "<<nIter<<"; deltaT 1;\n";
	ofs << "writeControl timeStep; writeInterval "<<nIter<<";\n";
	ofs << "purgeWrite 0; writeFormat ascii; writePrecision 8; writeCompression off;\n";
	ofs << "timeFormat general; timePrecision 8; runTimeModifiable true; libs (\"libatmosphericModels.so\");\n";
	ofs.close();
}

void Simulation::schemes(std::string directory) {
	std::string directoryS = directory+"/system/fvSchemes";
	std::ofstream ofs;
	ofs.open (directoryS.c_str());
	ofs << "FoamFile{version 2.0; format ascii; class dictionary; location \"system\"; object fvSchemes;}\n";
	ofs << "ddtSchemes{default steadyState;} gradSchemes{default Gauss linear;}\n";
	ofs << "divSchemes{default none; div(phi,U) bounded Gauss linearUpwindV Gauss linear;\n";
	ofs << "div(phi,k) bounded Gauss linearUpwind Gauss linear; div(phi,epsilon) bounded Gauss linearUpwind Gauss linear;\n";
	ofs << "div(phi,omega) bounded Gauss linearUpwind Gauss linear; div((nuEff*dev2(T(grad(U))))) Gauss linear;\n";
	ofs << "div(((nu+nut)*dev2(T(grad(U))))) Gauss linear;}\n";
	ofs << "laplacianSchemes{default Gauss linear corrected;} interpolationSchemes{default linear;}\n";
	ofs << "snGradSchemes{default corrected;} fluxRequired{default no; p;}";
	ofs.close();
}

void Simulation::a_schemes(std::string directory) {
	std::string directoryS = directory+"/system/fvSchemes";
	std::ofstream ofs;
	ofs.open (directoryS.c_str());
	ofs << "FoamFile{version 2.0; format ascii; class dictionary; location \"system\"; object fvSchemes;}\n";
	ofs << "ddtSchemes{default steadyState;} gradSchemes{default Gauss linear;}\n";
	ofs << "divSchemes{default none; div(-phi,U_a) bounded Gauss linearUpwindV Gauss linear;\n";
	ofs << "div(-phi,k_a) bounded Gauss linearUpwind Gauss linear; div(-phi,epsilon_a) bounded Gauss linearUpwind Gauss linear;\n";
	ofs << "div(-phi,omega_a) bounded Gauss linearUpwind Gauss linear;\n";
	ofs << "div(((nu+nut)*dev2(T(grad(U_a))))) Gauss linear; div((omega_a*dev2(T(grad(U))))) Gauss linear;\n";
	ofs << "div((((k*k_a)|omega)*symm(grad(U)))) Gauss linear; div((((k_a*k)|omega)*dev2(T(grad(U))))) Gauss linear;\n";
	ofs << "div((omega_a*symm(grad(U)))) Gauss linear;}\n";
	ofs << "laplacianSchemes{default Gauss linear corrected;} interpolationSchemes{default linear;}\n";
	ofs << "snGradSchemes{default corrected;} fluxRequired{default no; p; p_a;}";
	ofs.close();
}

void Simulation::solvers(std::string directory) {
	std::string directoryS = directory+"/system/fvSolution";
	std::ofstream ofs;
	ofs.open (directoryS.c_str());
	ofs << "FoamFile{version 2.0; format ascii; class dictionary; location \"system\"; object fvSolution;}\n";
	ofs << "solvers{\"(p|Phi)\"{solver GAMG; tolerance 1e-07; relTol 0.01; smoother GaussSeidel; nPreSweeps 0; nPostSweeps 2;\n";
	ofs << "cacheAgglomeration on; agglomerator faceAreaPair; nCellsInCoarsestLevel 200; mergeLevels 1;}\n";
	ofs << "\"(U|k|epsilon|omega)\"{solver smoothSolver; smoother GaussSeidel; tolerance 1e-08; relTol 0.01; nSweeps 2;}}\n";
	ofs << "SIMPLE{nNonOrthogonalCorrectors 0; residualControl{p 1e-5; \"(U|k|epsilon|omega)\" 1e-6;}}\n";
	ofs << "relaxationFactors{fields{p 0.3;} equations{\"(U|k)\" 0.5; \"(epsilon|omega)\" 0.3;}} cache{grad(U);}\n";
	ofs << "potentialFlow{nNonOrthogonalCorrectors 10;}";
	ofs.close();
}

void Simulation::a_solvers(int turb, std::string directory) {
	std::string directoryS = directory+"/system/fvSolution";
	std::ofstream ofs;
	ofs.open (directoryS.c_str());
	ofs << "FoamFile{version 2.0; format ascii; class dictionary; location \"system\"; object fvSolution;}\n";
	ofs << "solvers{p_a{solver GAMG; tolerance 1e-07; relTol 0.01; smoother GaussSeidel; nPreSweeps 0; nPostSweeps 2;\n";
	ofs << "cacheAgglomeration on; agglomerator faceAreaPair; nCellsInCoarsestLevel 200; mergeLevels 1;}\n";
	ofs << "\"(U_a|k_a|epsilon_a|omega_a)\"{solver smoothSolver; smoother GaussSeidel; tolerance 1e-08; relTol 0.01; nSweeps 2;}}\n";
	if (turb == 0) 
	{
		ofs << "SIMPLE{nNonOrthogonalCorrectors 0; residualControl{p_a 5e-5; \"(U_a|k_a|epsilon_a|omega_a)\" 5e-5;} p_aRefCell 0; p_aRefValue 0;}\n";
	}
	else 
	{
		ofs << "SIMPLE{nNonOrthogonalCorrectors 0; residualControl{p_a 1e-5; \"(U_a|k_a|epsilon_a|omega_a)\" 1e-6;} p_aRefCell 0; p_aRefValue 0;}\n";
	}
	ofs << "relaxationFactors{fields{p_a 0.3;} equations{\"(U_a|k_a|epsilon_a|omega_a)\" 0.5;}} cache{grad(U_a);}";
	ofs.close();
}

void Simulation::topoSet(std::vector< std::vector<double> > positions, double d, double h, double delta, int dim, std::string directory) {
	int rows =  positions.size();
	std::string directoryS = directory+"/system/topoSetDict";
	std::ofstream ofs;
	ofs.open (directoryS.c_str());
	ofs << "FoamFile{version 2.0; format ascii; class dictionary; location \"system\"; object topoSetDict;}\n";
	ofs << "actions (";
	if (dim == 2) 
	{
		for (int i=0; i<rows; i++)
		{
			ofs << "{name disk"<<i+1<<"; type cellSet; action new; source boxToCell; ";
			ofs << "sourceInfo{box ("<<positions[i][0]-d/2<<" "<<positions[i][1]-d/20<<" 0) ("<<positions[i][0]+d/2<<" "<<positions[i][1]+d/20<<" 1);}}\n";
			ofs << "{name disk"<<i+1<<"; type cellZoneSet; action new; source setToCellZone; sourceInfo{set disk"<<i+1<<";}}\n";
			ofs << "{name disk"<<i+1<<"plusX; type cellSet; action new; source boxToCell; ";
			ofs << "sourceInfo{box ("<<positions[i][0]-d/2+delta<<" "<<positions[i][1]-d/20<<" 0) ("<<positions[i][0]+d/2+delta<<" "<<positions[i][1]+d/20<<" 1);}}\n";
			ofs << "{name disk"<<i+1<<"plusX; type cellZoneSet; action new; source setToCellZone; sourceInfo{set disk"<<i+1<<"plusX;}}\n";
			ofs << "{name disk"<<i+1<<"minusX; type cellSet; action new; source boxToCell; ";
			ofs << "sourceInfo{box ("<<positions[i][0]-d/2-delta<<" "<<positions[i][1]-d/20<<" 0) ("<<positions[i][0]+d/2-delta<<" "<<positions[i][1]+d/20<<" 1);}}\n";
			ofs << "{name disk"<<i+1<<"minusX; type cellZoneSet; action new; source setToCellZone; sourceInfo{set disk"<<i+1<<"minusX;}}\n";
			ofs << "{name disk"<<i+1<<"plusY; type cellSet; action new; source boxToCell; ";
			ofs << "sourceInfo{box ("<<positions[i][0]-d/2<<" "<<positions[i][1]-d/20+delta<<" 0) ("<<positions[i][0]+d/2<<" "<<positions[i][1]+d/20+delta<<" 1);}}\n";
			ofs << "{name disk"<<i+1<<"plusY; type cellZoneSet; action new; source setToCellZone; sourceInfo{set disk"<<i+1<<"plusY;}}\n";
			ofs << "{name disk"<<i+1<<"minusY; type cellSet; action new; source boxToCell; ";
			ofs << "sourceInfo{box ("<<positions[i][0]-d/2<<" "<<positions[i][1]-d/20-delta<<" 0) ("<<positions[i][0]+d/2<<" "<<positions[i][1]+d/20-delta<<" 1);}}\n";
			ofs << "{name disk"<<i+1<<"minusY; type cellZoneSet; action new; source setToCellZone; sourceInfo{set disk"<<i+1<<"minusY;}}\n";
		}
	}
	else if (dim == 3)
	{
		for (int i=0; i<rows; i++)
		{
			ofs << "{name disk"<<i+1<<"; type cellSet; action new; source cylinderToCell; ";
			ofs << "sourceInfo{p1 ("<<positions[i][0]<<" "<<positions[i][1]-d/40<<" "<<h<<"); p2 ("<<positions[i][0]<<" "<<positions[i][1]+d/40<<" "<<h<<"); radius "<<d/2<<";}}\n";
			ofs << "{name disk"<<i+1<<"; type cellZoneSet; action new; source setToCellZone; sourceInfo{set disk"<<i+1<<";}}\n";
			ofs << "{name disk"<<i+1<<"plusX; type cellSet; action new; source cylinderToCell; ";
			ofs << "sourceInfo{p1 ("<<positions[i][0]+delta<<" "<<positions[i][1]-d/40<<" "<<h<<"); p2 ("<<positions[i][0]+delta<<" "<<positions[i][1]+d/40<<" "<<h<<"); radius "<<d/2<<";}}\n";
			ofs << "{name disk"<<i+1<<"plusX; type cellZoneSet; action new; source setToCellZone; sourceInfo{set disk"<<i+1<<"plusX;}}\n";
			ofs << "{name disk"<<i+1<<"minusX; type cellSet; action new; source cylinderToCell; ";
			ofs << "sourceInfo{p1 ("<<positions[i][0]-delta<<" "<<positions[i][1]-d/40<<" "<<h<<"); p2 ("<<positions[i][0]-delta<<" "<<positions[i][1]+d/40<<" "<<h<<"); radius "<<d/2<<";}}\n";
			ofs << "{name disk"<<i+1<<"minusX; type cellZoneSet; action new; source setToCellZone; sourceInfo{set disk"<<i+1<<"minusX;}}\n";
			ofs << "{name disk"<<i+1<<"plusY; type cellSet; action new; source cylinderToCell; ";
			ofs << "sourceInfo{p1 ("<<positions[i][0]<<" "<<positions[i][1]-d/40+delta<<" "<<h<<"); p2 ("<<positions[i][0]<<" "<<positions[i][1]+d/40+delta<<" "<<h<<"); radius "<<d/2<<";}}\n";
			ofs << "{name disk"<<i+1<<"plusY; type cellZoneSet; action new; source setToCellZone; sourceInfo{set disk"<<i+1<<"plusY;}}\n";
			ofs << "{name disk"<<i+1<<"minusY; type cellSet; action new; source cylinderToCell; ";
			ofs << "sourceInfo{p1 ("<<positions[i][0]<<" "<<positions[i][1]-d/40-delta<<" "<<h<<"); p2 ("<<positions[i][0]<<" "<<positions[i][1]+d/40-delta<<" "<<h<<"); radius "<<d/2<<";}}\n";
			ofs << "{name disk"<<i+1<<"minusY; type cellZoneSet; action new; source setToCellZone; sourceInfo{set disk"<<i+1<<"minusY;}}\n";
		}
	}
	ofs << ");";
	ofs.close();
}

void Simulation::topoSetComplexTerrain(double d, double h, double delta, std::string directory) {
	std::vector< std::vector<double> > turbines;
	std::vector<double> temp(3);
	std::string line;
	std::istringstream iss;
	std::string directoryS = directory+"/turbinesRotated.txt";
	std::ifstream inputs (directoryS.c_str());
	std::getline(inputs, line);
	int rows = 0;
	while (std::getline(inputs, line))
	{
		iss.str(line);
		iss >> temp[0];
		iss >> temp[1];
		iss >> temp[2];
		turbines.push_back(temp);
		rows = rows+1;
		iss.clear();
	}
	inputs.close();
	std::vector<double> slopeY(rows);
	directoryS = directory+"/terrainRotatedStreamWiseSlope.txt";
	inputs.open(directoryS.c_str());
	std::getline(inputs, line);
	for (int i=0; i<rows; i++)
	{
		std::getline(inputs, line);
		iss.str(line);
		iss >> slopeY[i];
		iss.clear();
	}
	inputs.close();
	std::vector<double> slopeX(rows);
	directoryS = directory+"/terrainRotatedCrossWiseSlope.txt";
	inputs.open(directoryS.c_str());
	std::getline(inputs, line);
	for (int i=0; i<rows; i++)
	{
		std::getline(inputs, line);
		iss.str(line);
		iss >> slopeX[i];
		iss.clear();
	}
	inputs.close();
	directoryS = directory+"/system/topoSetDict";
	std::ofstream ofs;
	ofs.open (directoryS.c_str());
	ofs << "FoamFile{version 2.0; format ascii; class dictionary; location \"system\"; object topoSetDict;}\n";
	ofs << "actions (";
	for (int i=0; i<rows; i++)
	{
		ofs << "{name disk"<<i+1<<"; type cellSet; action new; source cylinderToCell; ";
		ofs << "sourceInfo{p1 ("<<turbines[i][0]<<" "<<turbines[i][1]-d/40<<" "<<turbines[i][2]+h<<"); p2 ("<<turbines[i][0]<<" "<<turbines[i][1]+d/40<<" "<<turbines[i][2]+h<<"); radius "<<d/2<<";}}\n";
		ofs << "{name disk"<<i+1<<"; type cellZoneSet; action new; source setToCellZone; sourceInfo{set disk"<<i+1<<";}}\n";
		ofs << "{name disk"<<i+1<<"plusX; type cellSet; action new; source cylinderToCell; ";
		ofs << "sourceInfo{p1 ("<<turbines[i][0]+delta<<" "<<turbines[i][1]-d/40<<" "<<turbines[i][2]+slopeX[i]*delta+h<<"); p2 ("<<turbines[i][0]+delta<<" "<<turbines[i][1]+d/40<<" "<<turbines[i][2]+slopeX[i]*delta+h<<"); radius "<<d/2<<";}}\n";
		ofs << "{name disk"<<i+1<<"plusX; type cellZoneSet; action new; source setToCellZone; sourceInfo{set disk"<<i+1<<"plusX;}}\n";
		ofs << "{name disk"<<i+1<<"minusX; type cellSet; action new; source cylinderToCell; ";
		ofs << "sourceInfo{p1 ("<<turbines[i][0]-delta<<" "<<turbines[i][1]-d/40<<" "<<turbines[i][2]-slopeX[i]*delta+h<<"); p2 ("<<turbines[i][0]-delta<<" "<<turbines[i][1]+d/40<<" "<<turbines[i][2]-slopeX[i]*delta+h<<"); radius "<<d/2<<";}}\n";
		ofs << "{name disk"<<i+1<<"minusX; type cellZoneSet; action new; source setToCellZone; sourceInfo{set disk"<<i+1<<"minusX;}}\n";
		ofs << "{name disk"<<i+1<<"plusY; type cellSet; action new; source cylinderToCell; ";
		ofs << "sourceInfo{p1 ("<<turbines[i][0]<<" "<<turbines[i][1]-d/40+delta<<" "<<turbines[i][2]+slopeY[i]*delta+h<<"); p2 ("<<turbines[i][0]<<" "<<turbines[i][1]+d/40+delta<<" "<<turbines[i][2]+slopeY[i]*delta+h<<"); radius "<<d/2<<";}}\n";
		ofs << "{name disk"<<i+1<<"plusY; type cellZoneSet; action new; source setToCellZone; sourceInfo{set disk"<<i+1<<"plusY;}}\n";
		ofs << "{name disk"<<i+1<<"minusY; type cellSet; action new; source cylinderToCell; ";
		ofs << "sourceInfo{p1 ("<<turbines[i][0]<<" "<<turbines[i][1]-d/40-delta<<" "<<turbines[i][2]-slopeY[i]*delta+h<<"); p2 ("<<turbines[i][0]<<" "<<turbines[i][1]+d/40-delta<<" "<<turbines[i][2]-slopeY[i]*delta+h<<"); radius "<<d/2<<";}}\n";
		ofs << "{name disk"<<i+1<<"minusY; type cellZoneSet; action new; source setToCellZone; sourceInfo{set disk"<<i+1<<"minusY;}}\n";
	}
	ofs << ");";
	ofs.close();
}

void Simulation::a_topoSet(std::vector< std::vector<double> > positions, double lx, double ly, double lz, double d, double h, double delta, int dim, std::string directory) {
	int rows =  positions.size();
	std::string directoryS = directory+"/system/topoSetDict_a";
	std::ofstream ofs;
	ofs.open (directoryS.c_str());
	ofs << "FoamFile{version 2.0; format ascii; class dictionary; location \"system\"; object topoSetDict;}\n";
	ofs << "actions (";
	if (dim == 2) 
	{
		for (int i=0; i<rows; i++)
		{
			ofs << "{name disk"<<i+1<<"; type cellSet; action new; source boxToCell; ";
			ofs << "sourceInfo{box ("<<positions[i][0]-d/2<<" "<<positions[i][1]-d/20<<" 0) ("<<positions[i][0]+d/2<<" "<<positions[i][1]+d/20<<" 1);}}\n";
			ofs << "{name disk"<<i+1<<"; type cellZoneSet; action new; source setToCellZone; sourceInfo{set disk"<<i+1<<";}}\n";
			ofs << "{name disk"<<i+1<<"plusX; type cellSet; action new; source boxToCell; ";
			ofs << "sourceInfo{box ("<<positions[i][0]+d/2<<" "<<positions[i][1]-d/20<<" 0) ("<<positions[i][0]+d/2+delta<<" "<<positions[i][1]+d/20<<" 1);}}\n";
			ofs << "{name disk"<<i+1<<"plusX; type cellZoneSet; action new; source setToCellZone; sourceInfo{set disk"<<i+1<<"plusX;}}\n";
			ofs << "{name disk"<<i+1<<"plusXX; type cellSet; action new; source boxToCell; ";
			ofs << "sourceInfo{box ("<<positions[i][0]-d/2<<" "<<positions[i][1]-d/20<<" 0) ("<<positions[i][0]-d/2+delta<<" "<<positions[i][1]+d/20<<" 1);}}\n";
			ofs << "{name disk"<<i+1<<"plusXX; type cellZoneSet; action new; source setToCellZone; sourceInfo{set disk"<<i+1<<"plusXX;}}\n";
			ofs << "{name disk"<<i+1<<"minusX; type cellSet; action new; source boxToCell; ";
			ofs << "sourceInfo{box ("<<positions[i][0]+d/2-delta<<" "<<positions[i][1]-d/20<<" 0) ("<<positions[i][0]+d/2<<" "<<positions[i][1]+d/20<<" 1);}}\n";
			ofs << "{name disk"<<i+1<<"minusX; type cellZoneSet; action new; source setToCellZone; sourceInfo{set disk"<<i+1<<"minusX;}}\n";
			ofs << "{name disk"<<i+1<<"minusXX; type cellSet; action new; source boxToCell; ";
			ofs << "sourceInfo{box ("<<positions[i][0]-d/2-delta<<" "<<positions[i][1]-d/20<<" 0) ("<<positions[i][0]-d/2<<" "<<positions[i][1]+d/20<<" 1);}}\n";
			ofs << "{name disk"<<i+1<<"minusXX; type cellZoneSet; action new; source setToCellZone; sourceInfo{set disk"<<i+1<<"minusXX;}}\n";
			ofs << "{name disk"<<i+1<<"plusY; type cellSet; action new; source boxToCell; ";
			ofs << "sourceInfo{box ("<<positions[i][0]-d/2<<" "<<positions[i][1]+d/20<<" 0) ("<<positions[i][0]+d/2<<" "<<positions[i][1]+d/20+delta<<" 1);}}\n";
			ofs << "{name disk"<<i+1<<"plusY; type cellZoneSet; action new; source setToCellZone; sourceInfo{set disk"<<i+1<<"plusY;}}\n";
			ofs << "{name disk"<<i+1<<"plusYY; type cellSet; action new; source boxToCell; ";
			ofs << "sourceInfo{box ("<<positions[i][0]-d/2<<" "<<positions[i][1]-d/20<<" 0) ("<<positions[i][0]+d/2<<" "<<positions[i][1]-d/20+delta<<" 1);}}\n";
			ofs << "{name disk"<<i+1<<"plusYY; type cellZoneSet; action new; source setToCellZone; sourceInfo{set disk"<<i+1<<"plusYY;}}\n";
			ofs << "{name disk"<<i+1<<"minusY; type cellSet; action new; source boxToCell; ";
			ofs << "sourceInfo{box ("<<positions[i][0]-d/2<<" "<<positions[i][1]+d/20-delta<<" 0) ("<<positions[i][0]+d/2<<" "<<positions[i][1]+d/20<<" 1);}}\n";
			ofs << "{name disk"<<i+1<<"minusY; type cellZoneSet; action new; source setToCellZone; sourceInfo{set disk"<<i+1<<"minusY;}}\n";
			ofs << "{name disk"<<i+1<<"minusYY; type cellSet; action new; source boxToCell; ";
			ofs << "sourceInfo{box ("<<positions[i][0]-d/2<<" "<<positions[i][1]-d/20-delta<<" 0) ("<<positions[i][0]+d/2<<" "<<positions[i][1]-d/20<<" 1);}}\n";
			ofs << "{name disk"<<i+1<<"minusYY; type cellZoneSet; action new; source setToCellZone; sourceInfo{set disk"<<i+1<<"minusYY;}}\n";
		}
	}
	else if (dim == 3) 
	{
		for (int i=0; i<rows; i++)
		{
			ofs << "{name disk"<<i+1<<"; type cellSet; action new; source cylinderToCell; ";
			ofs << "sourceInfo{p1 ("<<positions[i][0]<<" "<<positions[i][1]-d/40<<" "<<h<<"); p2 ("<<positions[i][0]<<" "<<positions[i][1]+d/40<<" "<<h<<"); radius "<<d/2<<";}}\n";
			ofs << "{name disk"<<i+1<<"; type cellZoneSet; action new; source setToCellZone; sourceInfo{set disk"<<i+1<<";}}\n";
			ofs << "{name disk"<<i+1<<"plusX; type cellSet; action new; source cylinderToCell; ";
			ofs << "sourceInfo{p1 ("<<positions[i][0]+delta<<" "<<positions[i][1]-d/40<<" "<<h<<"); p2 ("<<positions[i][0]+delta<<" "<<positions[i][1]+d/40<<" "<<h<<"); radius "<<d/2<<";}}\n";
			ofs << "{name disk"<<i+1<<"plusX; type cellSet; action delete; source cellToCell; sourceInfo{set disk"<<i+1<<";}}\n";
			ofs << "{name disk"<<i+1<<"plusX; type cellZoneSet; action new; source setToCellZone; sourceInfo{set disk"<<i+1<<"plusX;}}\n";
			ofs << "{name disk"<<i+1<<"plusXX; type cellSet; action new; source cellToCell; sourceInfo{set disk"<<i+1<<";}}\n";
			ofs << "{name disk"<<i+1<<"plusXX; type cellSet; action delete; source cylinderToCell; ";
			ofs << "sourceInfo{p1 ("<<positions[i][0]+delta<<" "<<positions[i][1]-d/40<<" "<<h<<"); p2 ("<<positions[i][0]+delta<<" "<<positions[i][1]+d/40<<" "<<h<<"); radius "<<d/2<<";}}\n";
			ofs << "{name disk"<<i+1<<"plusXX; type cellZoneSet; action new; source setToCellZone; sourceInfo{set disk"<<i+1<<"plusXX;}}\n";
			ofs << "{name disk"<<i+1<<"minusX; type cellSet; action new; source cellToCell; sourceInfo{set disk"<<i+1<<";}}\n";
			ofs << "{name disk"<<i+1<<"minusX; type cellSet; action delete; source cylinderToCell; ";
			ofs << "sourceInfo{p1 ("<<positions[i][0]-delta<<" "<<positions[i][1]-d/40<<" "<<h<<"); p2 ("<<positions[i][0]-delta<<" "<<positions[i][1]+d/40<<" "<<h<<"); radius "<<d/2<<";}}\n";
			ofs << "{name disk"<<i+1<<"minusX; type cellZoneSet; action new; source setToCellZone; sourceInfo{set disk"<<i+1<<"minusX;}}\n";
			ofs << "{name disk"<<i+1<<"minusXX; type cellSet; action new; source cylinderToCell; ";
			ofs << "sourceInfo{p1 ("<<positions[i][0]-delta<<" "<<positions[i][1]-d/40<<" "<<h<<"); p2 ("<<positions[i][0]-delta<<" "<<positions[i][1]+d/40<<" "<<h<<"); radius "<<d/2<<";}}\n";
			ofs << "{name disk"<<i+1<<"minusXX; type cellSet; action delete; source cellToCell; sourceInfo{set disk"<<i+1<<";}}\n";
			ofs << "{name disk"<<i+1<<"minusXX; type cellZoneSet; action new; source setToCellZone; sourceInfo{set disk"<<i+1<<"minusXX;}}\n";
			ofs << "{name disk"<<i+1<<"plusY; type cellSet; action new; source cylinderToCell; ";
			ofs << "sourceInfo{p1 ("<<positions[i][0]<<" "<<positions[i][1]-d/40+delta<<" "<<h<<"); p2 ("<<positions[i][0]<<" "<<positions[i][1]+d/40+delta<<" "<<h<<"); radius "<<d/2<<";}}\n";
			ofs << "{name disk"<<i+1<<"plusY; type cellSet; action delete; source cellToCell; sourceInfo{set disk"<<i+1<<";}}\n";
			ofs << "{name disk"<<i+1<<"plusY; type cellZoneSet; action new; source setToCellZone; sourceInfo{set disk"<<i+1<<"plusY;}}\n";
			ofs << "{name disk"<<i+1<<"plusYY; type cellSet; action new; source cellToCell; sourceInfo{set disk"<<i+1<<";}}\n";
			ofs << "{name disk"<<i+1<<"plusYY; type cellSet; action delete; source cylinderToCell; ";
			ofs << "sourceInfo{p1 ("<<positions[i][0]<<" "<<positions[i][1]-d/40+delta<<" "<<h<<"); p2 ("<<positions[i][0]<<" "<<positions[i][1]+d/40+delta<<" "<<h<<"); radius "<<d/2<<";}}\n";
			ofs << "{name disk"<<i+1<<"plusYY; type cellZoneSet; action new; source setToCellZone; sourceInfo{set disk"<<i+1<<"plusYY;}}\n";
			ofs << "{name disk"<<i+1<<"minusY; type cellSet; action new; source cellToCell; sourceInfo{set disk"<<i+1<<";}}\n";
			ofs << "{name disk"<<i+1<<"minusY; type cellSet; action delete; source cylinderToCell; ";
			ofs << "sourceInfo{p1 ("<<positions[i][0]<<" "<<positions[i][1]-d/40-delta<<" "<<h<<"); p2 ("<<positions[i][0]<<" "<<positions[i][1]+d/40-delta<<" "<<h<<"); radius "<<d/2<<";}}\n";
			ofs << "{name disk"<<i+1<<"minusY; type cellZoneSet; action new; source setToCellZone; sourceInfo{set disk"<<i+1<<"minusY;}}\n";
			ofs << "{name disk"<<i+1<<"minusYY; type cellSet; action new; source cylinderToCell; ";
			ofs << "sourceInfo{p1 ("<<positions[i][0]<<" "<<positions[i][1]-d/40-delta<<" "<<h<<"); p2 ("<<positions[i][0]<<" "<<positions[i][1]+d/40-delta<<" "<<h<<"); radius "<<d/2<<";}}\n";
			ofs << "{name disk"<<i+1<<"minusYY; type cellSet; action delete; source cellToCell; sourceInfo{set disk"<<i+1<<";}}\n";
			ofs << "{name disk"<<i+1<<"minusYY; type cellZoneSet; action new; source setToCellZone; sourceInfo{set disk"<<i+1<<"minusYY;}}\n";
		}
		ofs << "{name toKeep; type cellSet; action new; source boxToCell; ";
		ofs << "sourceInfo{box ("<<-lx/2<<" "<<-ly/2<<" "<<d/2*0.005<<") ("<<lx/2<<" "<<ly/2<<" "<<lz<<");}}\n"; // to delete the first two cell layers at the wall
	}
	ofs << ");";
	ofs.close();
}

void Simulation::a_topoSetCompleTerrain(double lx, double ly, double lz, double d, double h, double delta, int dim, std::string directory) {
	std::vector< std::vector<double> > turbines;
	std::vector<double> temp(3);
	std::string line;
	std::istringstream iss;
	std::string directoryS = directory+"/turbinesRotated.txt";
	std::ifstream inputs (directoryS.c_str());
	std::getline(inputs, line);
	int rows = 0;
	while (std::getline(inputs, line))
	{
		iss.str(line);
		iss >> temp[0];
		iss >> temp[1];
		iss >> temp[2];
		turbines.push_back(temp);
		rows = rows+1;
		iss.clear();
	}
	inputs.close();
	std::vector<double> slopeY(rows);
	directoryS = directory+"/terrainRotatedStreamWiseSlope.txt";
	inputs.open(directoryS.c_str());
	std::getline(inputs, line);
	for (int i=0; i<rows; i++)
	{
		std::getline(inputs, line);
		iss.str(line);
		iss >> slopeY[i];
		iss.clear();
	}
	inputs.close();
	std::vector<double> slopeX(rows);
	directoryS = directory+"/terrainRotatedCrossWiseSlope.txt";
	inputs.open(directoryS.c_str());
	std::getline(inputs, line);
	for (int i=0; i<rows; i++)
	{
		std::getline(inputs, line);
		iss.str(line);
		iss >> slopeX[i];
		iss.clear();
	}
	inputs.close();
	directoryS = directory+"/system/topoSetDict_a";
	std::ofstream ofs;
	ofs.open (directoryS.c_str());
	ofs << "FoamFile{version 2.0; format ascii; class dictionary; location \"system\"; object topoSetDict;}\n";
	ofs << "actions (";
	for (int i=0; i<rows; i++)
	{
		ofs << "{name disk"<<i+1<<"; type cellSet; action new; source cylinderToCell; ";
		ofs << "sourceInfo{p1 ("<<turbines[i][0]<<" "<<turbines[i][1]-d/40<<" "<<turbines[i][2]+h<<"); p2 ("<<turbines[i][0]<<" "<<turbines[i][1]+d/40<<" "<<turbines[i][2]+h<<"); radius "<<d/2<<";}}\n";
		ofs << "{name disk"<<i+1<<"; type cellZoneSet; action new; source setToCellZone; sourceInfo{set disk"<<i+1<<";}}\n";
		ofs << "{name disk"<<i+1<<"plusX; type cellSet; action new; source cylinderToCell; ";
		ofs << "sourceInfo{p1 ("<<turbines[i][0]+delta<<" "<<turbines[i][1]-d/40<<" "<<turbines[i][2]+slopeX[i]*delta+h<<"); p2 ("<<turbines[i][0]+delta<<" "<<turbines[i][1]+d/40<<" "<<turbines[i][2]+slopeX[i]*delta+h<<"); radius "<<d/2<<";}}\n";
		ofs << "{name disk"<<i+1<<"plusX; type cellSet; action delete; source cellToCell; sourceInfo{set disk"<<i+1<<";}}\n";
		ofs << "{name disk"<<i+1<<"plusX; type cellZoneSet; action new; source setToCellZone; sourceInfo{set disk"<<i+1<<"plusX;}}\n";
		ofs << "{name disk"<<i+1<<"plusXX; type cellSet; action new; source cellToCell; sourceInfo{set disk"<<i+1<<";}}\n";
		ofs << "{name disk"<<i+1<<"plusXX; type cellSet; action delete; source cylinderToCell; ";
		ofs << "sourceInfo{p1 ("<<turbines[i][0]+delta<<" "<<turbines[i][1]-d/40<<" "<<turbines[i][2]+slopeX[i]*delta+h<<"); p2 ("<<turbines[i][0]+delta<<" "<<turbines[i][1]+d/40<<" "<<turbines[i][2]+slopeX[i]*delta+h<<"); radius "<<d/2<<";}}\n";
		ofs << "{name disk"<<i+1<<"plusXX; type cellZoneSet; action new; source setToCellZone; sourceInfo{set disk"<<i+1<<"plusXX;}}\n";
		ofs << "{name disk"<<i+1<<"minusX; type cellSet; action new; source cellToCell; sourceInfo{set disk"<<i+1<<";}}\n";
		ofs << "{name disk"<<i+1<<"minusX; type cellSet; action delete; source cylinderToCell; ";
		ofs << "sourceInfo{p1 ("<<turbines[i][0]-delta<<" "<<turbines[i][1]-d/40<<" "<<turbines[i][2]-slopeX[i]*delta+h<<"); p2 ("<<turbines[i][0]-delta<<" "<<turbines[i][1]+d/40<<" "<<turbines[i][2]-slopeX[i]*delta+h<<"); radius "<<d/2<<";}}\n";
		ofs << "{name disk"<<i+1<<"minusX; type cellZoneSet; action new; source setToCellZone; sourceInfo{set disk"<<i+1<<"minusX;}}\n";
		ofs << "{name disk"<<i+1<<"minusXX; type cellSet; action new; source cylinderToCell; ";
		ofs << "sourceInfo{p1 ("<<turbines[i][0]-delta<<" "<<turbines[i][1]-d/40<<" "<<turbines[i][2]-slopeX[i]*delta+h<<"); p2 ("<<turbines[i][0]-delta<<" "<<turbines[i][1]+d/40<<" "<<turbines[i][2]-slopeX[i]*delta+h<<"); radius "<<d/2<<";}}\n";
		ofs << "{name disk"<<i+1<<"minusXX; type cellSet; action delete; source cellToCell; sourceInfo{set disk"<<i+1<<";}}\n";
		ofs << "{name disk"<<i+1<<"minusXX; type cellZoneSet; action new; source setToCellZone; sourceInfo{set disk"<<i+1<<"minusXX;}}\n";
		ofs << "{name disk"<<i+1<<"plusY; type cellSet; action new; source cylinderToCell; ";
		ofs << "sourceInfo{p1 ("<<turbines[i][0]<<" "<<turbines[i][1]-d/40+delta<<" "<<turbines[i][2]+slopeY[i]*delta+h<<"); p2 ("<<turbines[i][0]<<" "<<turbines[i][1]+d/40+delta<<" "<<turbines[i][2]+slopeY[i]*delta+h<<"); radius "<<d/2<<";}}\n";
		ofs << "{name disk"<<i+1<<"plusY; type cellSet; action delete; source cellToCell; sourceInfo{set disk"<<i+1<<";}}\n";
		ofs << "{name disk"<<i+1<<"plusY; type cellZoneSet; action new; source setToCellZone; sourceInfo{set disk"<<i+1<<"plusY;}}\n";
		ofs << "{name disk"<<i+1<<"plusYY; type cellSet; action new; source cellToCell; sourceInfo{set disk"<<i+1<<";}}\n";
		ofs << "{name disk"<<i+1<<"plusYY; type cellSet; action delete; source cylinderToCell; ";
		ofs << "sourceInfo{p1 ("<<turbines[i][0]<<" "<<turbines[i][1]-d/40+delta<<" "<<turbines[i][2]+slopeY[i]*delta+h<<"); p2 ("<<turbines[i][0]<<" "<<turbines[i][1]+d/40+delta<<" "<<turbines[i][2]+slopeY[i]*delta+h<<"); radius "<<d/2<<";}}\n";
		ofs << "{name disk"<<i+1<<"plusYY; type cellZoneSet; action new; source setToCellZone; sourceInfo{set disk"<<i+1<<"plusYY;}}\n";
		ofs << "{name disk"<<i+1<<"minusY; type cellSet; action new; source cellToCell; sourceInfo{set disk"<<i+1<<";}}\n";
		ofs << "{name disk"<<i+1<<"minusY; type cellSet; action delete; source cylinderToCell; ";
		ofs << "sourceInfo{p1 ("<<turbines[i][0]<<" "<<turbines[i][1]-d/40-delta<<" "<<turbines[i][2]-slopeY[i]*delta+h<<"); p2 ("<<turbines[i][0]<<" "<<turbines[i][1]+d/40-delta<<" "<<turbines[i][2]-slopeY[i]*delta+h<<"); radius "<<d/2<<";}}\n";
		ofs << "{name disk"<<i+1<<"minusY; type cellZoneSet; action new; source setToCellZone; sourceInfo{set disk"<<i+1<<"minusY;}}\n";
		ofs << "{name disk"<<i+1<<"minusYY; type cellSet; action new; source cylinderToCell; ";
		ofs << "sourceInfo{p1 ("<<turbines[i][0]<<" "<<turbines[i][1]-d/40-delta<<" "<<turbines[i][2]-slopeY[i]*delta+h<<"); p2 ("<<turbines[i][0]<<" "<<turbines[i][1]+d/40-delta<<" "<<turbines[i][2]-slopeY[i]*delta+h<<"); radius "<<d/2<<";}}\n";
		ofs << "{name disk"<<i+1<<"minusYY; type cellSet; action delete; source cellToCell; sourceInfo{set disk"<<i+1<<";}}\n";
		ofs << "{name disk"<<i+1<<"minusYY; type cellZoneSet; action new; source setToCellZone; sourceInfo{set disk"<<i+1<<"minusYY;}}\n";
	}
	ofs << "{name toKeep; type cellSet; action new; source boxToCell; "; // to delete the first cell layer at the wall and at the inlet (inlet name will be replaced with terrain name)
	ofs << "sourceInfo{box ("<<-lx/2<<" "<<-ly/2<<" "<<0<<") ("<<lx/2<<" "<<ly/2<<" "<<lz<<");}}\n";
	ofs << "{name terrainFaces; type faceSet; action new; source patchToFace; sourceInfo{name terrain;}}\n";
	ofs << "{name terrainCells; type cellSet; action new; source faceToCell; sourceInfo{set terrainFaces; option owner;}}\n";
	ofs << "{name toKeep; type cellSet; action delete; source cellToCell; sourceInfo{set terrainCells;}}\n";
	ofs << "{name inletFaces; type faceSet; action new; source patchToFace; sourceInfo{name inlet;}}\n";
	ofs << "{name inletCells; type cellSet; action new; source faceToCell; sourceInfo{set inletFaces; option owner;}}\n";
	ofs << "{name toKeep; type cellSet; action delete; source cellToCell; sourceInfo{set inletCells;}}\n";
	ofs << ");";
	ofs.close();
}

void Simulation::options(int nwt, double T, std::string directory) {
	std::string directoryS = directory+"/constant/fvOptions";
	std::ofstream ofs;
	ofs.open (directoryS.c_str());
	ofs << "FoamFile{version 2.0; format ascii; class dictionary; location \"constant\"; object fvOptions;}\n";
	for (int i=0; i<nwt; i++)
	{
		ofs << "momentumSource"<<i+1<<"{type vectorSemiImplicitSource;\n";
		ofs << "vectorSemiImplicitSourceCoeffs{selectionMode cellZone; cellZone disk"<<i+1<<"; volumeMode absolute; injectionRateSuSp{U ((0 "<<T<<" 0) 0);}}}\n";
	}
	ofs.close();
}

void Simulation::a_options(int nwt, double T, std::string directory) {
	std::string directoryS = directory+"/constant/fvOptions";
	std::ofstream ofs;
	ofs.open (directoryS.c_str());
	ofs << "FoamFile{version 2.0; format ascii; class dictionary; location \"constant\"; object fvOptions;}\n";
	for (int i=0; i<nwt; i++)
	{
		ofs << "momentumSource"<<i+1<<"{type vectorSemiImplicitSource;\n";
		ofs << "vectorSemiImplicitSourceCoeffs{selectionMode cellZone; cellZone disk"<<i+1<<"; volumeMode absolute; injectionRateSuSp{U_a ((0 "<<T<<" 0) 0);}}}\n";
	}
	ofs.close();
}

void Simulation::run(int np, int turb, int dim, std::string directory) {
	std::cout << "Case " << directory << " (primal)\n";
	std::string command;
	runFoamCommand(directory,"blockMesh","",1,false);
	if (np>1) {
		std::ostringstream procSS;
		procSS << np;
		decompose(np,dim,directory);
		if (dim == 2) {
			command = "cp -r " + directory + "/0.primal/. " + directory + "/0/.";
			runCommand(command);
		}
		runFoamCommand(directory,"decomposePar","",1,false);
		if (dim == 3) {
			runFoamCommand(directory,"snappyHexMesh","-overwrite",np,false);
			cleanMesh(directory,np);
			addProcInterface(directory);
			command = "cd " + directory + "; ls -d processor* | xargs -I {} cp -r ./0.primal ./{}/0";
			runCommand(command);
		}
		runFoamCommand(directory,"renumberMesh","-overwrite",np,false);
		cleanMesh(directory,np);
		runFoamCommand(directory,"topoSet","",np,false);
		runFoamCommand(directory,"checkMesh","-allTopology -allGeometry",np,false);
		if (turb > 0) {
			if (dim == 3) {
				runFoamCommand(directory,"potentialFoam","",np,false);
			}
			runFoamCommand(directory,"simpleFoamS","",np,true);
		}
		else {
			runFoamCommand(directory,"icoFoamS","",np,true);
		}
		if (dim == 3) {
			runFoamCommand(directory,"reconstructParMesh","-constant",1,false);
		}
		runFoamCommand(directory,"reconstructPar","",1,false);
		command = "cd " + directory + "; ls -d processor* | xargs -I {} rm -rf ./{}";
		runCommand(command);
		cleanMesh(directory,1);
		runFoamCommand(directory,"topoSet","",1,false);
	}
	else {
		command = "cp -r " + directory + "/0.primal/. " + directory + "/0/.";
		runCommand(command);
		if (dim == 3) {
			runFoamCommand(directory,"snappyHexMesh","-overwrite",np,false);
			cleanMesh(directory,np);
		}
		runFoamCommand(directory,"topoSet","",np,false);
		runFoamCommand(directory,"checkMesh","",np,false);
		if (turb > 0) {
			if (dim == 3) {
				runFoamCommand(directory,"potentialFoam","",np,false);
			}
			runFoamCommand(directory,"simpleFoamS","",np,true);
		}
		else {
			runFoamCommand(directory,"icoFoamS","",np,true);
		}
	}
}

void Simulation::a_run(int np, int turb, int dim, std::string directory) {
	std::cout << "Case " << directory << " (adjoint)\n";
	std::string command;
	if (np>1) {
		command = "cp -r " + directory + "/0.adjoint/. " + directory + "/0/.";
		runCommand(command);
		command = "cp -r " + directory + "/constant.backup/. " + directory + "/constant/.";
		runCommand(command);
		decompose(np,dim,directory);
		runFoamCommand(directory,"decomposePar","",1,false);
		std::ostringstream procSS;
		procSS << np;
		if (dim == 3) {
			runFoamCommand(directory,"checkMesh","-allTopology -allGeometry",np,false);
			runFoamCommand(directory,"topoSet","-dict ../system/topoSetDict_a",np,false);
			runFoamCommand(directory,"subsetMesh","toKeep -patch terrain -overwrite",np,false);
			cleanMesh(directory,np);
			runFoamCommand(directory,"checkMesh","-allTopology -allGeometry",np,false);
		}
		runFoamCommand(directory,"renumberMesh","-overwrite",np,false);
		cleanMesh(directory,np);
		runFoamCommand(directory,"topoSet","-dict ../system/topoSetDict_a",np,false);
		if (turb == 1) {
			runFoamCommand(directory,"adjointSimpleFoamT","",np,true);
		}
		else if (turb == 2) {
			runFoamCommand(directory,"adjointSimpleFoam","",np,true);
		}
		else {
			runFoamCommand(directory,"adjointIcoFoam","",np,true);
		}
		runCommand(command);
		if (dim == 3) {
			runFoamCommand(directory,"reconstructParMesh","-constant",1,false);
		}
		runFoamCommand(directory,"reconstructPar","",1,false);
		command = "cd " + directory + "; ls -d processor* | xargs -I {} rm -rf ./{}";
		runCommand(command);
		cleanMesh(directory,1);
		runFoamCommand(directory,"topoSet","",1,false);
		command = "cd " + directory + "/constant/polyMesh; mv -f ./sets ./sets-init";
		runCommand(command);
		runFoamCommand(directory,"topoSet","-dict ./system/topoSetDict_a",1,false);
	}
	else {
		command = "cp -r " + directory + "/0.adjoint/. " + directory + "/0/.";
		runCommand(command);
		command = "cp -r " + directory + "/constant.backup/. " + directory + "/constant/.";
		runCommand(command);
		if (dim == 3) {
			runFoamCommand(directory,"checkMesh","",np,false);
			runFoamCommand(directory,"topoSet","-dict ./system/topoSetDict_a",np,false);
			runFoamCommand(directory,"subsetMesh","toKeep -patch terrain -overwrite",np,false);
			cleanMesh(directory,np);
			runFoamCommand(directory,"topoSet","",np,false);
			command = "cd " + directory + "/constant/polyMesh; mv -f ./sets ./sets-init";
			runCommand(command);
			command = directory + "/constant/polyMesh/cellZones";
			remove(command.c_str());
			runFoamCommand(directory,"checkMesh","",np,false);
		}
		runFoamCommand(directory,"topoSet","-dict ./system/topoSetDict_a",np,false);
		if (turb == 1) {
			runFoamCommand(directory,"adjointSimpleFoamT","",np,true);
		}
		else if (turb == 2) {
			runFoamCommand(directory,"adjointSimpleFoam","",np,true);
		}
		else {
			runFoamCommand(directory,"adjointIcoFoam","",np,true);
		}
	}
}

void Simulation::decompose(int np, int dim, std::string directory) {
	int nx;
	int ny;
	int nz;
	int rest1 = 100000;
	int rest2 = 100000;
	int trial1 = 0;
	int trial2 = 0;
	if (dim == 2)
	{
		while (rest1 != 0)
		{
			nx = floor(pow(np,0.5)) - trial1;
			trial1 = trial1 + 1;
			rest1 = remainder(np,nx);
		}
		ny = np/nx;
	}
	else if (dim == 3)
	{
		while (rest1 != 0)
		{
			nx = floor(pow(np,0.33333)) - trial1;
			trial1 = trial1 + 1;
			rest1 = remainder(np,nx);
		}
		while (rest2 != 0)
		{
			ny = floor(pow(np/nx,0.5)) - trial2;
			trial2 = trial2 + 1;
			rest2 = remainder(np/nx,ny);
		}
		nz = np/(nx*ny);
	}
	std::string directoryS = directory+"/system/decomposeParDict";
	std::ofstream ofs;
	ofs.open (directoryS.c_str());
	ofs << "FoamFile{version 2.0; format ascii; class dictionary; location \"system\"; object decomposeParDict;}\n";
	ofs << "numberOfSubdomains "<<np<<"; method simple; simpleCoeffs{n (";
	if (dim == 2)
	{
		ofs <<nx<<" "<<ny<<" 1";
	}
	else if (dim == 3)
	{
		ofs <<nx<<" "<<ny<<" "<<nz;
	}
	ofs << "); delta 0.001;}\n";
	ofs.close();
}

void Simulation::cleanMesh(std::string directory, int np) {
	if (np>1) {
		std::ostringstream procSS;
		for (int i=0; i<np; i++) {
			procSS.str("");
			procSS << i;
			std::string str = directory + "/processor" + procSS.str() + "/constant/polyMesh/surfaceIndex";
			remove(str.c_str());
			str = directory + "/processor" + procSS.str() + "/constant/polyMesh/refinementHistory";
			remove(str.c_str());
			str = directory + "/processor" + procSS.str() + "/constant/polyMesh/pointZones";
			remove(str.c_str());
			str = directory + "/processor" + procSS.str() + "/constant/polyMesh/pointLevel";
			remove(str.c_str());
			str = directory + "/processor" + procSS.str() + "/constant/polyMesh/level0Edge";
			remove(str.c_str());
			str = directory + "/processor" + procSS.str() + "/constant/polyMesh/faceZones";
			remove(str.c_str());
			str = directory + "/processor" + procSS.str() + "/constant/polyMesh/cellZones";
			remove(str.c_str());
			str = directory + "/processor" + procSS.str() + "/constant/polyMesh/cellLevel";
			remove(str.c_str());
			str = "cd " + directory + "/processor" + procSS.str() + "/constant/polyMesh; rm -rf ./sets";
			runCommand(str);
			str = "cd " + directory + "/processor" + procSS.str() + "/constant/polyMesh; rm -rf ./sets-init";
			runCommand(str);
		}

	}
	else {
		std::string str = directory + "/constant/polyMesh/surfaceIndex";
		remove(str.c_str());
		str = directory + "/constant/polyMesh/refinementHistory";
		remove(str.c_str());
		str = directory + "/constant/polyMesh/pointZones";
		remove(str.c_str());
		str = directory + "/constant/polyMesh/pointLevel";
		remove(str.c_str());
		str = directory + "/constant/polyMesh/level0Edge";
		remove(str.c_str());
		str = directory + "/constant/polyMesh/faceZones";
		remove(str.c_str());
		str = directory + "/constant/polyMesh/cellZones";
		remove(str.c_str());
		str = directory + "/constant/polyMesh/cellLevel";
		remove(str.c_str());
		str = "cd " + directory + "/constant/polyMesh; rm -rf ./sets";
		runCommand(str);
		str = "cd " + directory + "/constant/polyMesh; rm -rf ./sets-init";
		runCommand(str);
	}
}

void Simulation::cleanPrimal(std::string directory) {
	std::string str = directory + "/PrimalResultFiles";
	mkdir(str.c_str(), S_IRWXU | S_IRWXG | S_IROTH | S_IXOTH);
	str = "cd " + directory + "; ls -d *Output | xargs -I {} mv {} ./PrimalResultFiles/.";
	runCommand(str);
	str = "cd " + directory + "; ls -d core* | xargs -I {} rm {}";
	runCommand(str);
}

void Simulation::cleanAdjoint(std::string directory) {
	std::string str = directory + "/AdjointResultFiles";
	mkdir(str.c_str(), S_IRWXU | S_IRWXG | S_IROTH | S_IXOTH);
	str = "cd " + directory + "; ls -d *Output | xargs -I {} mv {} ./AdjointResultFiles/.";
	runCommand(str);
	str = "cd " + directory + "; rm -r ./constant.backup";
	runCommand(str);
	str = "cd " + directory + "; ls -d core* | xargs -I {} rm {}";
	runCommand(str);
}

void Simulation::setAdjoint(int finalIter, std::string directory) {
	std::string str = "cd " + directory + "; mv ./0 ./PrimalResultFiles/.";
	runCommand(str);
	str = "cd " + directory + "; mv ./system ./PrimalResultFiles/.";
	runCommand(str);
	str = directory + "/constant/RASProperties";
	remove(str.c_str());
	str = directory + "/constant/fvOptions";
	remove(str.c_str());
	str = directory + "/constant/turbulenceProperties";
	remove(str.c_str());
	std::ostringstream oss;
	oss << "cd " << directory << "; mv ./" << finalIter << " ./0.adjoint";
	str = oss.str();
	runCommand(str);
	str = "cd " + directory + "/0.adjoint; rm -r ./uniform";
	runCommand(str);
	str = "cd " + directory + "/constant/polyMesh; rm -r ./cellZones ";
	runCommand(str);
	str = "cd " + directory + "/constant/polyMesh; mv ./sets ./sets-init";
	runCommand(str);
	str = "cd " + directory + "; rm -r ./dynamicCode";
	runCommand(str);
	str = "cd " + directory + "; mv ./constant ./constant.backup";
	runCommand(str);
}

void Simulation::addProcInterface(std::string directory) {
	std::vector<std::string> temp;
	std::string line;
	std::ifstream filein;
	std::ofstream fileout;
	std::vector<std::string> directoryS;

	directoryS.push_back(directory+"/0.primal/U");
	directoryS.push_back(directory+"/0.primal/p");
	directoryS.push_back(directory+"/0.primal/k");
	directoryS.push_back(directory+"/0.primal/omega");
	directoryS.push_back(directory+"/0.primal/epsilon");
	directoryS.push_back(directory+"/0.primal/nut");
	for (int i = 0; i < directoryS.size(); i++)
	{
		filein.open(directoryS[i].c_str());
		int n = 0;
		while (std::getline(filein, line))
		{
			temp.push_back(line);
			n = n+1;
		}
		filein.close();
		fileout.open(directoryS[i].c_str());
		for (int j = 0; j < n-1; j++)
		{
			fileout << temp[j] << std::endl;
		}
		fileout << "sideW{type zeroGradient;}" << std::endl;
		fileout << "\"proc.*\"{type processor;}}";
		fileout.close();
		temp.clear();
	}
	directoryS.clear();
}

void Simulation::runCommand(std::string command) {
	std::string commandToRun = "timeout 30m bash -c '" + command + "'";
	int checkSistem = 0;
	while (checkSistem < 10000) {
		if (system(NULL)) {
			system(commandToRun.c_str());
			checkSistem = 10000;
		}
		else {
			std::cout << "Waiting for system availability...";
			checkSistem = checkSistem + 1;
		}
	}
}

void Simulation::runFoamCommand(std::string directory, std::string command, std::string options, int np, bool userDefined) {
	std::string OFdir = "";
	std::string OFUdir = "";
	//std::string OFUdir = "/home/c/camon/enrico8/OpenFOAM/enrico8-3.0.1/platforms/linux64GccDPInt32Opt/bin/";
	//std::string OFUdir = "/home/c/camon/enrico8/OpenFOAM/enrico8-5.0/platforms/linux64GccDPInt32Opt/bin/";
	//std::string OFUdir = "/Carnegie/DGE/Homes/Users/eantonini/OpenFOAM/eantonini-6/platforms/linux64GccDPInt64Opt/bin/";
        //std::string MPIdir = "/opt/sharcnet/openmpi/2.1.1/gcc630-std/bin/";
	std::string MPIdir = "";
	int delayed = 0;
	int time;
	if (command.find("simpleFoamS") != std::string::npos) {
		//time = 40;
		time = 180;
		//time = 540;
		delayed = 1;
	}
	else if (command.find("icoFoamS") != std::string::npos) {
		//time = 40;
		time = 180;
		//time = 540;
	}
	else if (command.find("adjointSimpleFoamT") != std::string::npos) {
		//time = 30;
		time = 180;
		//time = 540;
	}
	else if (command.find("adjointSimpleFoam") != std::string::npos) {
		//time = 30;
		time = 180;
		//time = 540;
	}
	else if (command.find("adjointIcoFoam") != std::string::npos) {
		//time = 30;
		time = 180;
		//time = 540;
	}
	else {
		//time = 10;
		time = 120;
		//time = 240;
	}
	std::string delayedCommand = "sed -i 's/fields{p 0.3;} equations{\"(U|k)\" 0.5; \"(epsilon|omega)\" 0.3;}/fields{p 0.4;} equations{\"(U|k)\" 0.7; \"(epsilon|omega)\" 0.5;}/g' " + directory +"/system/fvSolution";
	std::string finalCommand;
	if (np>1) {
		std::ostringstream procSS;
		procSS << np;
		if (userDefined) {
			finalCommand = "cd " + directory + "; " + MPIdir + "mpirun --mca mca_component_show_load_errors 0 --mca mpi_warn_on_fork 0 -np " + procSS.str() + " " + OFUdir + command + " " + options + " -parallel > " + command + "Output";
		}
		else {
			finalCommand = "cd " + directory + "; " + MPIdir + "mpirun --mca mca_component_show_load_errors 0 --mca mpi_warn_on_fork 0 -np " + procSS.str() + " " + OFdir + command + " " + options + " -parallel > " + command + "Output";
		}
	}
	else {
		if (userDefined) {
			finalCommand = "cd " + directory + "; " + OFUdir + command + " " + options + " > " + command + "Output";
		}
		else {
			finalCommand = "cd " + directory + "; " + OFdir + command + " " + options + " > " + command + "Output";
		}
	}
	std::ostringstream oss;
	oss << "timeout " << time << "m bash -c '" << finalCommand << "'"; // Timeout in case the program freezes
	std::string commandToRun;
	if (delayed == 0) {
		commandToRun = oss.str();
	}
	else {
		commandToRun = "(( sleep 20m; " + delayedCommand + ") & " + oss.str() + " )";
	}
	int success = 0;
	while (success < 10) {
		if (success == 0) {
			std::cout << "Running " << command << " ... ";
		}
		else {
			std::cout << "\nRerunning " << command << " ... ";
		}
		int checkSistem = 0;
		while (checkSistem < 10000) {
			if (system(NULL)) {
				system(commandToRun.c_str());
				checkSistem = 10000;
			}
			else {
				std::cout << "Waiting for system availability..." << "\n";
				checkSistem = checkSistem + 1;
			}
		}
		std::string commandOutputDir = directory + "/" + command + "Output";  // Check if command succeeded
		std::ifstream ifs (commandOutputDir.c_str());
		std::string line1; std::getline(ifs, line1);
		std::string line2; std::getline(ifs, line2);
		std::string line3; std::getline(ifs, line3);
		std::string line4; std::getline(ifs, line4);
		std::string temp;
		while (std::getline(ifs, temp)) {
			line1 = line2;
			line2 = line3;
			line3 = line4;
			line4 = temp;
		}
		std::istringstream tempSS;
		tempSS.str(line1);
		if (line1.find("End") != std::string::npos) {success = 10;}
		tempSS.clear();
		tempSS.str(line2);
		if (line2.find("End") != std::string::npos) {success = 10;}
		tempSS.clear();
		tempSS.str(line3);
		if (line3.find("End") != std::string::npos) {success = 10;}
		tempSS.clear();
		tempSS.str(line4);
		if (line4.find("End") != std::string::npos) {success = 10;}
		tempSS.clear();
		ifs.close();
		success = success + 1;
	}
	std::cout << "done" << "\n";
}
