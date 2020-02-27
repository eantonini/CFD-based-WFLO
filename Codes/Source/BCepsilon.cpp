#include "BCepsilon.h"
#include <fstream>
#include <string>

BCepsilon::BCepsilon() {
}

void BCepsilon::twodLayout(double epsilon, std::string directory) {
	std::string directoryS = directory+"/0.primal/epsilon";
	std::ofstream ofs;
	ofs.open (directoryS.c_str());
	ofs << "FoamFile{version 2.0; format ascii; class volScalarField; location \"0\"; object epsilon;}\n";
	ofs << "dimensions [0 2 -3 0 0 0 0]; internalField uniform "<<epsilon<<";\n";
	ofs << "boundaryField{\n";
	ofs << "inlet{type fixedValue; value uniform "<<epsilon<<";}\n";
	ofs << "outlet{type zeroGradient;}\n";
	ofs << "upAndDown{type empty;}\n";
	ofs << "sideE{type zeroGradient;}\n";
	ofs << "sideW{type zeroGradient;}}\n";
	ofs.close();
}

void BCepsilon::a_twodLayout(std::string directory) {
	std::string directoryS = directory+"/0.adjoint/epsilon_a";
	std::ofstream ofs;
	ofs.open (directoryS.c_str());
	ofs << "FoamFile{version 2.0; format ascii; class volScalarField; location \"0\"; object epsilon_a;}\n";
	ofs << "dimensions [0 0 0 0 0 0 0]; internalField uniform 0;\n";
	ofs << "boundaryField{\n";
	ofs << "inlet{type fixedValue; value uniform 0;}\n";
	ofs << "outlet{type fixedValue; value uniform 0;}\n";
	ofs << "upAndDown{type empty;}\n";
	ofs << "sideE{type fixedValue; value uniform 0;}\n";
	ofs << "sideW{type fixedValue; value uniform 0;}}\n";
	ofs.close();
}

void BCepsilon::threedLayout(double epsilon, double epsilonTop, double Vtop, double lz, double z0, std::string directory) {
	std::string directoryS = directory+"/0.primal/epsilon";
	std::ofstream ofs;
	ofs.open (directoryS.c_str());
	ofs << "FoamFile{version 2.0; format ascii; class volScalarField; location \"0\"; object epsilon;}\n";
	ofs << "dimensions [0 2 -3 0 0 0 0]; internalField uniform "<<epsilon<<";\n";
	ofs << "boundaryField{\n";
	ofs << "inlet{type atmBoundaryLayerInletEpsilon; Uref "<<Vtop<<"; Zref	"<<lz<<"; zDir (0 0 1); flowDir (0 -1 0);\n";
	ofs << "z0 uniform "<<z0<<"; zGround uniform 0.0; Cmu 0.0333; kappa 0.41; value $internalField;}\n";
	ofs << "outlet{type zeroGradient;}\n";
	ofs << "terrain{type epsilonWallFunction; value uniform "<<epsilon<<";}\n";
	ofs << "top{type fixedValue; value uniform "<<epsilonTop<<";}\n";
	ofs << "sideE{type zeroGradient;}\n";
	ofs << "sideW{type zeroGradient;}}\n";
	ofs.close();
}

void BCepsilon::a_threedLayout(std::string directory) {
	std::string directoryS = directory+"/0.adjoint/epsilon_a";
	std::ofstream ofs;
	ofs.open (directoryS.c_str());
	ofs << "FoamFile{version 2.0; format ascii; class volScalarField; location \"0\"; object epsilon_a;}\n";
	ofs << "dimensions [0 0 0 0 0 0 0]; internalField uniform 0;\n";
	ofs << "boundaryField{\n";
	ofs << "inlet{type fixedValue; value uniform 0;}\n";
	ofs << "outlet{type fixedValue; value uniform 0;}\n";
	ofs << "terrain{type fixedValue; value uniform 0;}\n";
	ofs << "top{type fixedValue; value uniform 0;}\n";
	ofs << "sideE{type fixedValue; value uniform 0;}\n";
	ofs << "sideW{type fixedValue; value uniform 0;}}\n";
	ofs.close();
}
