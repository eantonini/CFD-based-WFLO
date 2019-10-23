#ifndef SRC_BCEPSILON_H_
#define SRC_BCEPSILON_H_
#include <string>

class BCepsilon {
public:
	BCepsilon();
	void twodLayout(double epsilon, std::string directory);
	void a_twodLayout(std::string directory);
	void threedLayout(double epsilon, double epsilonMax, double Vmax, double lz, double z0, std::string directory);
	void a_threedLayout(std::string directory);
};

#endif
