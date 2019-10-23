#ifndef SRC_BCVELOCITY_H_
#define SRC_BCVELOCITY_H_
#include <string>

class BCvelocity {
public:
	BCvelocity();
	void twodLayout(double V, std::string directory);
	void a_twodLayout(std::string directory);
	void threedLayout(double V, double Vmax, double lz, double z0, std::string directory);
	void a_threedLayout(std::string directory);
};

#endif
