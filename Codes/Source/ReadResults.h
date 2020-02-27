#ifndef SRC_READRESULTS_H_
#define SRC_READRESULTS_H_
#include <vector>
#include <string>

class ReadResults {
public:
	ReadResults();
	void finalIter(int& finalIter, int turb, std::string directory);
	void a_finalIter(int& finalIter, int turb, std::string directory);
	void getIter(int& finalIter, std::string directoryS);
	void rotorVelocities(std::vector<double>& rotorVel, int finalIter, std::string directory);
	void rotorVolumes(std::vector<double>& rotorVol, int finalIter, std::string directory);
	void dummyRotorVelocities(std::vector<double>& dummyVel, int finalIter, std::string directory);
	void dummyRotorVolumes(std::vector<double>& dummyVol, int finalIter, std::string directory);
	void adjointVelocities(std::vector<double>& adjointVel, int finalIter, std::string directory);
	void adjointVolumes(std::vector<double>& adjointVol, int finalIter, std::string directory);
	void averageVelocity(double& averageVel, std::string directoryDisk, std::string directoryU, std::string directoryVol);
	void totalVolume(double& totVol, std::string directoryDisk, std::string directoryVol);
};

#endif
