#ifndef SRC_SIMULATION_H_
#define SRC_SIMULATION_H_

#include <vector>
#include <string>

class Simulation {
public:
	Simulation();
	void control(int nIter, int turb, std::string directory);
	void a_control(int nIter, int turb, std::string directory);
	void schemes(std::string directory);
	void a_schemes(std::string directory);
	void solvers(std::string directory);
	void a_solvers(int turb, std::string directory);
	void options(int nwt, double T, std::string directory);
	void a_options(int nwt, double T, std::string directory);
	void topoSet(std::vector< std::vector<double> > positions, double d, double h, double delta, int dim, std::string directory);
	void topoSetComplexTerrain(double d, double h, double delta, std::string directory);
	void a_topoSet(std::vector< std::vector<double> > positions, double lx, double ly, double lz, double d, double h, double delta, int dim, std::string directory);
	void a_topoSetCompleTerrain(double lx, double ly, double lz, double d, double h, double delta, int dim, std::string directory);
	void run(int np, int turb, int dim, std::string directory);
	void a_run(int np, int turb, int dim, std::string directory);
	void decompose(int np, int dim, std::string directory);
	void cleanMesh(std::string directory, int np);
	void cleanPrimal(std::string directory);
	void cleanAdjoint(std::string directory);
	void setAdjoint(int finalIter, std::string directory);
	void addProcInterface(std::string directory);
	void runCommand(std::string command);
	void runFoamCommand(std::string directory, std::string command, std::string options, int np, bool userDefined);
};

#endif
