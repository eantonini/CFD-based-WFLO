#include "BCk.h"
#include <fstream>
#include <string>

BCk::BCk() {
}

void BCk::twodLayout(double k, std::string directory) {
	std::string directoryS = directory+"/0.primal/k";
	std::ofstream ofs;
	ofs.open (directoryS.c_str());
	ofs << "FoamFile{version 2.0; format ascii; class volScalarField; location \"0\"; object k;}\n";
	ofs << "dimensions [0 2 -2 0 0 0 0]; internalField uniform "<<k<<";\n";
	ofs << "boundaryField{\n";
	ofs << "inlet{type fixedValue; value uniform "<<k<<";}\n";
	ofs << "outlet{type zeroGradient;}\n";
	ofs << "upAndDown{type empty;}\n";
	ofs << "sideE{type zeroGradient;}\n";
	ofs << "sideW{type zeroGradient;}}\n";
	ofs.close();
}

void BCk::a_twodLayout(std::string directory) {
	std::string directoryS = directory+"/0.adjoint/k_a";
	std::ofstream ofs;
	ofs.open (directoryS.c_str());
	ofs << "FoamFile{version 2.0; format ascii; class volScalarField; location \"0\"; object k_a;}\n";
	ofs << "dimensions [0 0 0 0 0 0 0]; internalField uniform 0;\n";
	ofs << "boundaryField{\n";
	ofs << "inlet{type fixedValue; value uniform 0;}\n";
	ofs << "outlet{type fixedValue; value uniform 0;}\n";
	ofs << "upAndDown{type empty;}\n";
	ofs << "sideE{type fixedValue; value uniform 0;}\n";
	ofs << "sideW{type fixedValue; value uniform 0;}}\n";
	ofs.close();
}

void BCk::threedLayout(double k, std::string directory) {
	std::string directoryS = directory+"/0.primal/k";
	std::ofstream ofs;
	ofs.open (directoryS.c_str());
	ofs << "FoamFile{version 2.0; format ascii; class volScalarField; location \"0\"; object k;}\n";
	ofs << "dimensions [0 2 -2 0 0 0 0]; internalField uniform "<<k<<";\n";
	ofs << "boundaryField{\n";
	ofs << "inlet{type fixedValue; value uniform "<<k<<";}\n";
	ofs << "outlet{type zeroGradient;}\n";
	ofs << "terrain{type kqRWallFunction; value uniform "<<k<<";}\n";
	ofs << "top{type fixedValue; value uniform "<<k<<";}\n";
	ofs << "sideE{type zeroGradient;}\n";
	ofs << "sideW{type zeroGradient;}}\n";
	ofs.close();
}

void BCk::a_threedLayout(std::string directory) {
	std::string directoryS = directory+"/0.adjoint/k_a";
	std::ofstream ofs;
	ofs.open (directoryS.c_str());
	ofs << "FoamFile{version 2.0; format ascii; class volScalarField; location \"0\"; object k_a;}\n";
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
