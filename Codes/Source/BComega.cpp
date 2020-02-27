#include "BComega.h"
#include <fstream>
#include <string>

BComega::BComega() {
}

void BComega::twodLayout(double omega, std::string directory) {
	std::string directoryS = directory+"/0.primal/omega";
	std::ofstream ofs;
	ofs.open (directoryS.c_str());
	ofs << "FoamFile{version 2.0; format ascii; class volScalarField; location \"0\"; object omega;}\n";
	ofs << "dimensions [0 0 -1 0 0 0 0]; internalField uniform "<<omega<<";\n";
	ofs << "boundaryField{\n";
	ofs << "inlet{type fixedValue; value uniform "<<omega<<";}\n";
	ofs << "outlet{type zeroGradient;}\n";
	ofs << "upAndDown{type empty;}\n";
	ofs << "sideE{type zeroGradient;}\n";
	ofs << "sideW{type zeroGradient;}}\n";
	ofs.close();
}

void BComega::a_twodLayout(std::string directory) {
	std::string directoryS = directory+"/0.adjoint/omega_a";
	std::ofstream ofs;
	ofs.open (directoryS.c_str());
	ofs << "FoamFile{version 2.0; format ascii; class volScalarField; location \"0\"; object omega_a;}\n";
	ofs << "dimensions [0 2 -1 0 0 0 0]; internalField uniform 0;\n";
	ofs << "boundaryField{\n";
	ofs << "inlet{type fixedValue; value uniform 0;}\n";
	ofs << "outlet{type fixedValue; value uniform 0;}\n";
	ofs << "upAndDown{type empty;}\n";
	ofs << "sideE{type fixedValue; value uniform 0;}\n";
	ofs << "sideW{type fixedValue; value uniform 0;}}\n";
	ofs.close();
}

void BComega::threedLayout(double omega, double omegaTop, double Ustar, double z0, std::string directory) {
	std::string directoryS = directory+"/0.primal/omega";
	std::ofstream ofs;
	ofs.open (directoryS.c_str());
	ofs << "FoamFile{version 2.0; format ascii; class volScalarField; location \"0\"; object omega;}\n";
	ofs << "dimensions [0 0 -1 0 0 0 0]; internalField uniform "<<omega<<";\n";
	ofs << "boundaryField{\n";
	ofs << "inlet{type codedFixedValue; value uniform "<<omega<<"; redirectType atmBoundaryLayerInletOmega;\n";
	ofs << "code #{scalar Ustar = "<<Ustar<<"; scalar kappa = 0.41; scalar Cmu = 0.0333; scalar z0 = "<<z0<<"; vector zDir = vector(0,0,1);\n";
	ofs << "scalarField zComp = patch().Cf()&zDir; scalarField value = Ustar/(kappa*pow(Cmu,0.5)*(z0+zComp)); operator == (value); #};}\n";
	ofs << "outlet{type zeroGradient;}\n";
	ofs << "terrain{type omegaWallFunction; value uniform "<<omega<<";}\n";
	ofs << "top{type fixedValue; value uniform "<<omegaTop<<";}\n";
	ofs << "sideE{type zeroGradient;}\n";
	ofs << "sideW{type zeroGradient;}}\n";
	ofs.close();
}

void BComega::a_threedLayout(std::string directory) {
	std::string directoryS = directory+"/0.adjoint/omega_a";
	std::ofstream ofs;
	ofs.open (directoryS.c_str());
	ofs << "FoamFile{version 2.0; format ascii; class volScalarField; location \"0\"; object omega_a;}\n";
	ofs << "dimensions [0 2 -1 0 0 0 0]; internalField uniform 0;\n";
	ofs << "boundaryField{\n";
	ofs << "inlet{type fixedValue; value uniform 0;}\n";
	ofs << "outlet{type fixedValue; value uniform 0;}\n";
	ofs << "terrain{type fixedValue; value uniform 0;}\n";
	ofs << "top{type fixedValue; value uniform 0;}\n";
	ofs << "sideE{type fixedValue; value uniform 0;}\n";
	ofs << "sideW{type fixedValue; value uniform 0;}}\n";
	ofs.close();
}
