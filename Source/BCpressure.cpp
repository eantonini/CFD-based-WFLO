#include "BCpressure.h"
#include <fstream>
#include <string>

BCpressure::BCpressure() {
}

void BCpressure::twodLayout(std::string directory) {
	std::string directoryS = directory+"/0.primal/p";
	std::ofstream ofs;
	ofs.open (directoryS.c_str());
	ofs << "FoamFile{version 2.0; format ascii; class volScalarField; location \"0\"; object p;}\n";
	ofs << "dimensions [0 2 -2 0 0 0 0]; internalField uniform 0;\n";
	ofs << "boundaryField{\n";
	ofs << "inlet{type zeroGradient;}\n";
	ofs << "outlet{type fixedValue; value uniform 0;}\n";
	ofs << "upAndDown{type empty;}\n";
	ofs << "sideE{type zeroGradient;}\n";
	ofs << "sideW{type zeroGradient;}}\n";
	ofs.close();
}

void BCpressure::a_twodLayout(std::string directory) {
	std::string directoryS = directory+"/0.adjoint/p_a";
	std::ofstream ofs;
	ofs.open (directoryS.c_str());
	ofs << "FoamFile{version 2.0; format ascii; class volScalarField; location \"0\"; object p_a;}\n";
	ofs << "dimensions [0 2 -2 0 0 0 0]; internalField uniform 0;\n";
	ofs << "boundaryField{\n";
	ofs << "inlet{type zeroGradient;}\n";
	ofs << "outlet{type zeroGradient;}\n";
	ofs << "upAndDown{type empty;}\n";
	ofs << "sideE{type zeroGradient;}\n";
	ofs << "sideW{type zeroGradient;}}\n";
	ofs.close();
}

void BCpressure::threedLayout(std::string directory) {
	std::string directoryS = directory+"/0.primal/p";
	std::ofstream ofs;
	ofs.open (directoryS.c_str());
	ofs << "FoamFile{version 2.0; format ascii; class volScalarField; location \"0\"; object p;}\n";
	ofs << "dimensions [0 2 -2 0 0 0 0]; internalField uniform 0;\n";
	ofs << "boundaryField{\n";
	ofs << "inlet{type zeroGradient;}\n";
	ofs << "outlet{type fixedValue; value uniform 0;}\n";
	ofs << "terrain{type zeroGradient;}\n";
	ofs << "top{type zeroGradient;}\n";
	ofs << "sideE{type zeroGradient;}\n";
	ofs << "sideW{type zeroGradient;}}\n";
	ofs.close();
}

void BCpressure::a_threedLayout(std::string directory) {
	std::string directoryS = directory+"/0.adjoint/p_a";
	std::ofstream ofs;
	ofs.open (directoryS.c_str());
	ofs << "FoamFile{version 2.0; format ascii; class volScalarField; location \"0\"; object p_a;}\n";
	ofs << "dimensions [0 2 -2 0 0 0 0]; internalField uniform 0;\n";
	ofs << "boundaryField{\n";
	ofs << "inlet{type zeroGradient;}\n";
	ofs << "outlet{type zeroGradient;}\n";
	ofs << "terrain{type zeroGradient;}\n";
	ofs << "top{type zeroGradient;}\n";
	ofs << "sideE{type zeroGradient;}\n";
	ofs << "sideW{type zeroGradient;}}\n";
	ofs.close();
}
