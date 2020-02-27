#include "ReadInputs.h"
#include <sstream>
#include <fstream>
#include <vector>
#include <iostream>

ReadInputs::ReadInputs() {
}

void ReadInputs::setup(double& d, double& h, double& h0, double& z0, double& lx, double& ly, double& lz, double& delta, int& turb, int& dim, int& np, int& terrain, double& height, double& spread) {
	std::ifstream inputs("InputsPowerCalculation.txt");
	std::string line;
	std::getline(inputs, line);
	std::istringstream iss;
	iss.str(line);
	iss >> d;
	iss.clear();
	std::getline(inputs, line);
	iss.str(line);
	iss >> h;
	iss.clear();
	std::getline(inputs, line);
	iss.str(line);
	iss >> h0;
	iss.clear();
	std::getline(inputs, line);
	iss.str(line);
	iss >> z0;
	iss.clear();
	std::getline(inputs, line);
	iss.str(line);
	iss >> lx;
	iss.clear();
	std::getline(inputs, line);
	iss.str(line);
	iss >> ly;
	iss.clear();
	std::getline(inputs, line);
	iss.str(line);
	iss >> lz;
	iss.clear();
	std::getline(inputs, line);
	iss.str(line);
	iss >> delta;
	iss.clear();
	std::getline(inputs, line);
	iss.str(line);
	iss >> turb;
	iss.clear();
	std::getline(inputs, line);
	iss.str(line);
	iss >> dim;
	iss.clear();
	std::getline(inputs, line);
	iss.str(line);
	iss >> np;
	iss.clear();
	std::getline(inputs, line);
	iss.str(line);
	iss >> terrain;
	iss.clear();
	std::getline(inputs, line);
	iss.str(line);
	iss >> height;
	iss.clear();
	std::getline(inputs, line);
	iss.str(line);
	iss >> spread;
	iss.clear();
}

void ReadInputs::thrustCoeff(std::vector<double>& Ct) {
	std::ifstream inputs("InputsThrustCoefficient.txt");
	std::string line;
	std::getline(inputs, line);
	std::istringstream iss;
	int rows = Ct.size();
	for (int i=0; i<rows; i++)
	{
		std::getline(inputs, line);
		iss.str(line);
		iss >> Ct[i];
		iss >> Ct[i];
		iss.clear();
	}
}


void ReadInputs::positions(int& nt, std::vector< std::vector<double> >& positions) {
	std::ifstream inputs("InputsPositions.txt");
	std::string line;
	std::getline(inputs, line);
	std::istringstream iss;
	std::vector<double> temp(2);
	nt = 0;
	while (std::getline(inputs, line))
	{
		iss.str(line);
		iss >> temp[0];
		iss >> temp[1];
		positions.push_back(temp);
		nt = nt+1;
		iss.clear();
	}
}

void ReadInputs::windRose(int wsbinmax, int& nwd, std::vector<double>& windDir, std::vector< std::vector<double> >& windRose) {
	std::ifstream inputs("InputsWindRose.txt");
	std::string line;
	std::getline(inputs, line);
	std::istringstream iss(line);
	double tempD;
	nwd = 0;
	while (iss >> tempD)
	{
		windDir.push_back(tempD);
		nwd = nwd+1;
	}
	iss.clear();
	for (int i=0; i<wsbinmax; i++)
	{
		std::vector<double> temp;
		std::getline(inputs, line);
		iss.str(line);
		iss >> tempD;
		while (iss >> tempD)
		{
			temp.push_back(tempD);
		}
		iss.clear();
		windRose.push_back(temp);
	}
}
