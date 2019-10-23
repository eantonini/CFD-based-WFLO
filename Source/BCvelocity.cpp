#include "BCvelocity.h"
#include <fstream>
#include <string>

BCvelocity::BCvelocity() {
}

void BCvelocity::twodLayout(double V, std::string directory) {
	std::string directoryS = directory+"/0.primal/U";
	std::ofstream ofs;
	ofs.open (directoryS.c_str());
	ofs << "FoamFile{version 2.0; format ascii; class volVectorField; location \"0\"; object U;}\n";
	ofs << "dimensions [0 1 -1 0 0 0 0]; internalField uniform (0 "<<-V<<" 0);\n";
	ofs << "boundaryField{\n";
	ofs << "inlet{type fixedValue; value uniform (0 "<<-V<<" 0);}\n";
	ofs << "outlet{type zeroGradient;}\n";
	ofs << "upAndDown{type empty;}\n";
	ofs << "sideE{type zeroGradient;}\n";
	ofs << "sideW{type zeroGradient;}}\n";
	ofs.close();
}

void BCvelocity::a_twodLayout(std::string directory) {
	std::string directoryS = directory+"/0.adjoint/U_a";
	std::ofstream ofs;
	ofs.open (directoryS.c_str());
	ofs << "FoamFile{version 2.0; format ascii; class volVectorField; location \"0\"; object U_a;}\n";
	ofs << "dimensions [0 1 -1 0 0 0 0]; internalField uniform (0 0 0);\n";
	ofs << "boundaryField{\n";
	ofs << "inlet{type fixedValue; value uniform (0 0 0);}\n";
	ofs << "outlet{type fixedValue; value uniform (0 0 0);}\n";
	ofs << "upAndDown{type empty;}\n";
	ofs << "sideE{type fixedValue; value uniform (0 0 0);}\n";
	ofs << "sideW{type fixedValue; value uniform (0 0 0);}}\n";
	ofs.close();
}

void BCvelocity::threedLayout(double V, double Vtop, double lz, double z0, std::string directory) {
	std::string directoryS = directory+"/0.primal/U";
	std::ofstream ofs;
	ofs.open (directoryS.c_str());
	ofs << "FoamFile{version 2.0; format ascii; class volVectorField; location \"0\"; object U;}\n";
	ofs << "dimensions [0 1 -1 0 0 0 0]; internalField uniform (0 "<<-V<<" 0);\n";
	ofs << "boundaryField{\n";
	//ofs << "inlet{type atmBoundaryLayerInletVelocity; Uref "<<Vtop<<"; Zref "<<lz<<"; zDir (0 0 1); flowDir (0 -1 0);\n";
	//ofs << "z0 uniform "<<z0<<"; zGround uniform 0.0; Cmu 0.0333; kappa 0.41; value $internalField;}\n";
	ofs << "inlet{type codedFixedValue; value uniform (0 "<<-V<<" 0); redirectType atmBoundaryLayerUDFInletVelocity;\n";
	ofs << "code #{scalar Uref = "<<Vtop<<"; scalar Zref = "<<lz<<"; scalar kappa = 0.41; scalar z0 = "<<z0<<"; vector zDir = vector(0,0,1);\n";
	ofs << "scalar Ustar = Uref*kappa/log((Zref+z0)/z0); vectorField value(patch().size(), vector(0,0,0));\n";
	ofs << "forAll(patch().Cf(),i) {scalar zComp = patch().Cf()[i]&zDir; value[i] = vector(0,-Ustar/kappa*log((z0+zComp)/z0),0);} operator == (value); #};}\n";
	ofs << "outlet{type zeroGradient;}\n";
	ofs << "terrain{type fixedValue; value uniform (0 0 0);}\n";
	ofs << "top{type fixedValue; value uniform (0 "<<-Vtop<<" 0);}\n";
	ofs << "sideE{type zeroGradient;}\n";
	ofs << "sideW{type zeroGradient;}}\n";
	ofs.close();
}

void BCvelocity::a_threedLayout(std::string directory) {
	std::string directoryS = directory+"/0.adjoint/U_a";
	std::ofstream ofs;
	ofs.open (directoryS.c_str());
	ofs << "FoamFile{version 2.0; format ascii; class volVectorField; location \"0\"; object U_a;}\n";
	ofs << "dimensions [0 1 -1 0 0 0 0]; internalField uniform (0 0 0);\n";
	ofs << "boundaryField{\n";
	ofs << "inlet{type fixedValue; value uniform (0 0 0);}\n";
	ofs << "outlet{type fixedValue; value uniform (0 0 0);}\n";
	ofs << "terrain{type fixedValue; value uniform (0 0 0);}\n";
	ofs << "top{type fixedValue; value uniform (0 0 0);}\n";
	ofs << "sideE{type fixedValue; value uniform (0 0 0);}\n";
	ofs << "sideW{type fixedValue; value uniform (0 0 0);}}\n";
	ofs.close();
}
