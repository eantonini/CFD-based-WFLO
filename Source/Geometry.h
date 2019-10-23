#ifndef SRC_GEOMETRY_H_
#define SRC_GEOMETRY_H_
#include <vector>
#include <string>

class Geometry {
public:
	Geometry();
	void twodLayout(double lx, double ly, double delta, std::string directory);
	void threedLayout(double lx, double ly, double lz, double d, double delta, std::string directory);
	void threedLayoutRefine(std::vector< std::vector<double> > positions, double d, double h, std::string directory);
	void gaussianHill(double hight, double spread, double domainWidth, double resolution);
	void remapTerrain(double windDir, std::string directory);
	void remapTurbines(std::vector< std::vector<double> > positions, std::string directory);
	void threedLayoutComplexTerrain(double lx, double ly, double lz, double d, double delta, std::string directory);
	void threedLayoutRefineComplexTerrain(double d, double h, std::string directory);
};

#endif
