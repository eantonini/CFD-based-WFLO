#include "DiffusionProperties.h"
#include <fstream>
#include <string>

DiffusionProperties::DiffusionProperties() {
}

void DiffusionProperties::turbulentKE(std::string directory) {
	std::string directoryS = directory+"/constant/turbulenceProperties";
	std::ofstream ofs;
	ofs.open (directoryS.c_str());
	ofs << "FoamFile{version 2.0; format ascii; class dictionary; location \"constant\"; object turbulenceProperties;}\n";
	ofs << "simulationType RAS;";
	ofs << "RAS{RASModel kEpsilon; turbulence on; printCoeffs off;}";
	ofs.close();
}

void DiffusionProperties::turbulentKO(std::string directory) {
	std::string directoryS = directory+"/constant/turbulenceProperties";
	std::ofstream ofs;
	ofs.open (directoryS.c_str());
	ofs << "FoamFile{version 2.0; format ascii; class dictionary; location \"constant\"; object turbulenceProperties;}\n";
	ofs << "simulationType RAS;";
	ofs << "RAS{RASModel kOmega; turbulence on; printCoeffs on;\n";
	//ofs << "kOmegaCoeffs{betaStar 0.0333; gamma 0.42; beta 0.0277; alphaK 0.45; alphaOmega 0.45;}}";
	ofs << "kOmegaCoeffs{betaStar 0.0333; gamma 1.16; beta 0.055; alphaK 0.54; alphaOmega 0.54;}}";
	ofs.close();
}

void DiffusionProperties::turbulentKOSST(std::string directory) {
	std::string directoryS = directory+"/constant/turbulenceProperties";
	std::ofstream ofs;
	ofs.open (directoryS.c_str());
	ofs << "FoamFile{version 2.0; format ascii; class dictionary; location \"constant\"; object turbulenceProperties;}\n";
	ofs << "simulationType RAS;";
	ofs << "RAS{RASModel kOmegaSST; turbulence on; printCoeffs off;}";
	ofs.close();
}

void DiffusionProperties::viscosity(double nu, std::string directory) {
	std::string directoryS = directory+"/constant/transportProperties";
	std::ofstream ofs;
	ofs.open (directoryS.c_str());
	ofs << "FoamFile{version 2.0; format ascii; class dictionary; location \"constant\"; object transportProperties;}\n";
	ofs << "transportModel Newtonian; nu nu [0 2 -1 0 0 0 0] "<<nu<<";";
	ofs.close();
}
