#ifndef SRC_READINPUTS_H_
#define SRC_READINPUTS_H_
#include <vector>

class ReadInputs {
public:
	ReadInputs();
	void setup(double& d, double& h, double& h0, double& z0, double& lx, double& ly, double& lz, double& delta, int& turb, int& dim, int& np, int& terrain, double& height, double& spread);
	void thrustCoeff(std::vector<double>& Ct);
	void positions(int& nt, std::vector< std::vector<double> >& positions);
	void windRose(int wsbinmax, int& nwd, std::vector<double>& windDir, std::vector< std::vector<double> >& windRose);
};

#endif
