#ifndef SRC_BCPRESSURE_H_
#define SRC_BCPRESSURE_H_
#include <string>

class BCpressure {
public:
	BCpressure();
	void twodLayout(std::string directory);
	void a_twodLayout(std::string directory);
	void threedLayout(std::string directory);
	void a_threedLayout(std::string directory);
};

#endif
