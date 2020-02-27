#ifndef SRC_DIFFUSIONPROPERTIES_H_
#define SRC_DIFFUSIONPROPERTIES_H_
#include <string>

class DiffusionProperties {
public:
	DiffusionProperties();
	void turbulentKE(std::string directory);
	void turbulentKO(std::string directory);
	void turbulentKOSST(std::string directory);
	void viscosity(double nu, std::string directory);
};

#endif
