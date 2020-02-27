#include "BCnut.h"
#include <fstream>
#include <string>

BCnut::BCnut() {
}

void BCnut::twodLayout(std::string directory) {
	std::string directoryS = directory+"/0.primal/nut";
	std::ofstream ofs;
	ofs.open (directoryS.c_str());
	ofs << "FoamFile{version 2.0; format ascii; class volScalarField; location \"0\"; object nut;}\n";
	ofs << "dimensions [0 2 -1 0 0 0 0]; internalField uniform 0;\n";
	ofs << "boundaryField{\n";
	ofs << "inlet{type calculated; value uniform 0;}\n";
	ofs << "outlet{type calculated; value uniform 0;}\n";
	ofs << "upAndDown{type empty;}\n";
	ofs << "sideE{type calculated; value uniform 0;}\n";
	ofs << "sideW{type calculated; value uniform 0;}}\n";
	ofs.close();
}

void BCnut::threedLayout(double z0, std::string directory) {
	std::string directoryS = directory+"/0.primal/nut";
	std::ofstream ofs;
	ofs.open (directoryS.c_str());
	ofs << "FoamFile{version 2.0; format ascii; class volScalarField; location \"0\"; object nut;}\n";
	ofs << "dimensions [0 2 -1 0 0 0 0]; internalField uniform 0;\n";
	ofs << "boundaryField{\n";
	ofs << "inlet{type calculated; value uniform 0;}\n";
	ofs << "outlet{type calculated; value uniform 0;}\n";
	ofs << "terrain{type nutkAtmRoughWallFunction; z0 uniform "<<z0<<"; value uniform 0.0;}\n";
	ofs << "top{type calculated; value uniform 0;}\n";
	ofs << "sideE{type calculated; value uniform 0;}\n";
	ofs << "sideW{type calculated; value uniform 0;}}\n";
	ofs.close();
}
