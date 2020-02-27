#ifndef SRC_BCNUT_H_
#define SRC_BCNUT_H_
#include <string>

class BCnut {
public:
	BCnut();
	void twodLayout(std::string directory);
	void threedLayout(double z0, std::string directory);
};

#endif
