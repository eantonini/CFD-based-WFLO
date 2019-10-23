#include <math.h>
#include <cmath>
#include <vector>
#include <sys/stat.h>
#include <sys/types.h>
#include <sstream>
#include <iostream>
#include <fstream>
#include <string>
#include <iomanip>

#include "nlopt.hpp"

#define pi 3.1416

int count = 0;
int nT;
double R;
double D;
std::ofstream resultFile;

double objFunction(const std::vector<double> &x, std::vector<double> &grad, void *obj_data)
{
	++count;
	std::ostringstream directorySS;
	std::string str;

	if (count>11) {
		directorySS << "rm -rf ./Iteration_" << count-10;
		str = directorySS.str();
		system(str.c_str());
		directorySS.str("");
	}

	double objValue = 0;
	
	directorySS << "./Iteration_" << count;
	std::string directoryS = directorySS.str();
	mkdir(directoryS.c_str(), S_IRWXU | S_IRWXG | S_IROTH | S_IXOTH);
	str = "cp -f ./PowerCalculation ./" + directoryS + "/";
	system(str.c_str());
	str = "cp -f ./GradientCalculation ./" + directoryS + "/";
	system(str.c_str());
	str = "cp -f ./RemapTerrain.py ./" + directoryS + "/";
	system(str.c_str());
	str = "cp -f ./RemapTurbines.py ./" + directoryS + "/";
	system(str.c_str());
	str = "cp -f ./InputsPowerCalculation.txt ./" + directoryS + "/";
	system(str.c_str());
	str = "cp -f ./InputsThrustCoefficient.txt ./" + directoryS + "/";
	system(str.c_str());
	str = "cp -f ./InputsWindRose.txt ./" + directoryS + "/";
	system(str.c_str());
	str = "cp -f ./terrain.txt ./" + directoryS + "/";
	system(str.c_str());
	str = "./" + directoryS +"/InputsPositions.txt";
	std::ofstream ofs;
	ofs.open (str.c_str());
	ofs << "x\ty\n";
	for (int i=0; i<nT; i++) {
		ofs << x[i*2+0] << "\t" << x[i*2+1] << "\n";
	}
	ofs.close();

	str = "cd ./" + directoryS +"; ./PowerCalculation";
	system(str.c_str());

	if (!grad.empty()) {
		str = "cd ./" + directoryS +"; ./GradientCalculation";
		system(str.c_str());
	}

	str = "./" + directoryS +"/Results-EnergyProduction.txt";
	std::ifstream ifs (str.c_str());
	std::string finalLine;
	std::string tempLine;
	while (std::getline(ifs, tempLine)) {
		finalLine = tempLine;
	}
	std::istringstream tempSS;
	tempSS.str(finalLine);
	tempSS.seekg(32);
	tempSS >> objValue;
	tempSS.clear();
	ifs.close();

	if (!grad.empty()) {
		str = "./" + directoryS +"/Results-Gradient.txt";
		ifs.open(str.c_str());
		std::string finalLines[nT*2+2];
		while (std::getline(ifs, tempLine)) {
			for (int i=0; i<nT*2+1; i++) {
				finalLines[i] = finalLines[i+1];
			}
			finalLines[nT*2+1] = tempLine;
		}
		if (!grad.empty()) {
			for (int i=0; i<nT; i++) {
				tempSS.str(finalLines[i]);
				tempSS >> grad[i*2+0];
				tempSS >> grad[i*2+1];
				tempSS.clear();
			}
		}
		ifs.close();
	}

	resultFile.open("Results.txt", std::ofstream::app);
	resultFile << std::right;
	resultFile << std::setw(6) << count;
	resultFile << std::fixed << std::setprecision(1);
	resultFile << std::setw(12) << objValue;
	for (int i=0; i<nT; i++) {
		resultFile << std::setw(12) << x[i*2+0] << std::setw(12) << x[i*2+1];
	}
	if (!grad.empty()) {
		for (int i=0; i<nT; i++) {
			resultFile << std::setw(12) << grad[i*2+0] << std::setw(12) << grad[i*2+1];
		}
	}
	else {
		for (int i=0; i<nT; i++) {
			resultFile << std::setw(12) << 0 << std::setw(12) << 0;
		}
	}
	resultFile << "\n";
	resultFile.close();

	return objValue;
}

void constraint1(unsigned m, double *result, unsigned n, const double* x, double* grad, void* constr_data) // Outer boundary constraint
{
	for (int i=0; i<nT; i++) {
		result[i] = pow(pow(x[i*2+0],2)+pow(x[i*2+1],2),2)-pow(R,4);
		for (int j=0; j<nT; j++) {
			if (grad) {
				if (j == i) {
					grad[i*nT*2+j*2+0] = 2*(pow(x[i*2+0],2)+pow(x[i*2+1],2))*2*x[i*2+0];
					grad[i*nT*2+j*2+1] = 2*(pow(x[i*2+0],2)+pow(x[i*2+1],2))*2*x[i*2+1];
				}
				else {
					grad[i*nT*2+j*2+0] = 0;
					grad[i*nT*2+j*2+1] = 0;
				}
			}
		}
	}
}

void constraint2(unsigned m, double *result, unsigned n, const double* x, double* grad, void* constr_data) // Turbine interspacing constraint
{
	for (int i=0; i<nT; i++) {
		for (int j=0; j<nT; j++) {
			if (j == i) {
				result[i*nT+j] = 0;
			}
			else {
				result[i*nT+j] = pow(1*D,2)-pow(x[i*2+0]-x[j*2+0],2)-pow(x[i*2+1]-x[j*2+1],2);
			}
			for (int k=0; k<nT; k++) {
				if (grad) {
					if (k == i) {
						grad[(i*nT+j)*nT*2+k*2+0] = -2*(x[i*2+0]-x[j*2+0]);
						grad[(i*nT+j)*nT*2+k*2+1] = -2*(x[i*2+1]-x[j*2+1]);
					}
					else if (k == j) {
						grad[(i*nT+j)*nT*2+k*2+0] = 2*(x[i*2+0]-x[j*2+0]);
						grad[(i*nT+j)*nT*2+k*2+1] = 2*(x[i*2+1]-x[j*2+1]);
					}
					else {
						grad[(i*nT+j)*nT*2+k*2+0] = 0;
						grad[(i*nT+j)*nT*2+k*2+1] = 0;
					}
				}
			}
		}
	}
}

int main() {

	std::istringstream iss;
	std::ostringstream directorySS;
	std::string str;
	std::string line;
	std::ifstream ifs;
	double temp;

	std::vector<double> x;

	int lastIter = 0;
	str = "./Results.txt";
	ifs.open(str.c_str());
	if (ifs.is_open()) {
		while (std::getline(ifs, line)) {
			iss.str(line);
			iss >> lastIter;
			iss.clear();
		}
	}
	ifs.close();
	
	count = count + lastIter;

	if (lastIter > 0) {
		directorySS << "./Iteration_" << lastIter+1 << "/InputsPositions.txt";
		str = directorySS.str();
		directorySS.str("");
	}
	else {
		str = "./InputsPositions.txt";
	}
	
	ifs.open(str.c_str());
	std::getline(ifs, line);
	nT = 0;
	while (std::getline(ifs, line)) {
		iss.str(line);
		iss >> temp;
		x.push_back(temp);
		iss >> temp;
		x.push_back(temp);
		nT = nT+1;
		iss.clear();
	}
	ifs.close();

	if (lastIter > 0) {
		directorySS << "rm -rf ./Iteration_" << lastIter+1;
		str = directorySS.str();
		directorySS.str("");
		system(str.c_str());
	}
	
	str = "./InputsPowerCalculation.txt";
	ifs.open(str.c_str());
	std::getline(ifs, line);
	iss.str(line);
	iss >> D;
	iss.clear();
	for (int i=0; i<4; i++) {
		std::getline(ifs, line);
	}
	iss.str(line);
	iss >> temp;
	R = temp/2 - 4*D; // Circular constraint
	ifs.close();

	nlopt::opt opt(nlopt::LD_SLSQP, nT*2);

	std::vector<double> lb(nT*2,-R);
	std::vector<double> ub(nT*2,R);
	opt.set_lower_bounds(lb);
	opt.set_upper_bounds(ub);

	opt.set_max_objective(objFunction, NULL);

	const std::vector<double> tol1(nT, 0);
	const std::vector<double> tol2(nT*nT, 0);
	opt.add_inequality_mconstraint(constraint1, NULL, tol1);
	opt.add_inequality_mconstraint(constraint2, NULL, tol2);

	opt.set_ftol_rel(0.01);

	double minf;
	nlopt::result result = opt.optimize(x, minf);

	resultFile.open("Results.txt", std::ofstream::app);
	if (minf < 0) {
		resultFile << "\nnlopt failed!";
	}
	else {
		resultFile << "\nFound minimum after " << count << " evaluations\n";
		resultFile << "at f(\n";
		for (int i=0; i<nT; i++) {
			resultFile << std::setw(12) << x[i*2+0] << ", " << std::setw(12) << x[i*2+1] << ", ";
		}
		resultFile << "\n) = " << minf;
	}
	resultFile.close();

	return 0;
}
