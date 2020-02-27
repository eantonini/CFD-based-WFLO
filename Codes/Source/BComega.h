#ifndef SRC_BCOMEGA_H_
#define SRC_BCOMEGA_H_
#include <string>

class BComega {
public:
	BComega();
	void twodLayout(double omega, std::string directory);
	void a_twodLayout(std::string directory);
	void threedLayout(double omega, double omegaMax, double Ustar, double z0, std::string directory);
	void a_threedLayout(std::string directory);
};

#endif
