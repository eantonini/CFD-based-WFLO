#include "ReadResults.h"
#include <sstream>
#include <fstream>
#include <vector>
#include <math.h>
#include <iostream>
#include <string>
#include <algorithm>

ReadResults::ReadResults() {
}

void ReadResults::finalIter(int& finalIter, int turb, std::string directory) {
	std::ostringstream directorySS;
	if (turb > 0) {
		directorySS << directory << "/simpleFoamSOutput";
	}
	else {
		directorySS << directory << "/icoFoamSOutput";
	}
	std::string directoryS = directorySS.str();
	getIter(finalIter,directoryS);
}

void ReadResults::a_finalIter(int& finalIter, int turb, std::string directory) {
	std::ostringstream directorySS;
	if (turb == 1) {
		directorySS << directory << "/adjointSimpleFoamTOutput";
	}
	else if (turb == 2) {
		directorySS << directory << "/adjointSimpleFoamOutput";
	}
	else {
		directorySS << directory << "/adjointIcoFoamOutput";
	}
	std::string directoryS = directorySS.str();
	getIter(finalIter,directoryS);
}

void ReadResults::getIter(int& finalIter, std::string directory) {
	std::ifstream ifs (directory.c_str());
	std::string line1; std::getline(ifs, line1);
	std::string line2; std::getline(ifs, line2);
	std::string line3; std::getline(ifs, line3);
	std::string line4; std::getline(ifs, line4);
	std::string line5; std::getline(ifs, line5);
	std::string temp;
	while (std::getline(ifs, temp)) {
		line1 = line2;
		line2 = line3;
		line3 = line4;
		line4 = line5;
		line5 = temp;
	}
	std::istringstream tempSS;
	std::string lineTemp = "F";
	if (line5[0] == lineTemp[0]) {  //Parallel
		tempSS.str(line1);
	}
	else {
		tempSS.str(line2);
	}
	tempSS.seekg(29);
	tempSS >> finalIter;
}

void ReadResults::rotorVelocities(std::vector<double>& rotorVel, int finalIter, std::string directory) {
	int rows =  rotorVel.size();
	std::ostringstream directorySS;
	directorySS << directory << "/" << finalIter << "/U";
	std::string directoryU = directorySS.str();
	directorySS.str("");
	directorySS << directory << "/" << finalIter << "/Vol";
	std::string directoryVol = directorySS.str();
	directorySS.str("");
	std::string directoryDisk;
	for (int i=0; i<rows; i++)
	{
		directorySS << directory << "/constant/polyMesh/sets/disk" << i+1;
		std::string directoryDisk = directorySS.str();
		directorySS.str("");
		averageVelocity(rotorVel[i],directoryDisk,directoryU,directoryVol);
	}
}

void ReadResults::rotorVolumes(std::vector<double>& rotorVol, int finalIter, std::string directory) {
	int rows =  rotorVol.size();
	std::ostringstream directorySS;
	directorySS << directory << "/" << finalIter << "/Vol";
	std::string directoryVol = directorySS.str();
	directorySS.str("");
	std::string directoryDisk;
	for (int i=0; i<rows; i++)
	{
		directorySS << directory << "/constant/polyMesh/sets/disk" << i+1;
		std::string directoryDisk = directorySS.str();
		directorySS.str("");
		totalVolume(rotorVol[i],directoryDisk,directoryVol);
	}
}

void ReadResults::dummyRotorVelocities(std::vector<double>& dummyVel, int finalIter, std::string directory) {
	int rows =  dummyVel.size()/4;
	std::ostringstream directorySS;
	directorySS << directory << "/" << finalIter << "/U";
	std::string directoryU = directorySS.str();
	directorySS.str("");
	directorySS << directory << "/" << finalIter << "/Vol";
	std::string directoryVol = directorySS.str();
	directorySS.str("");
	std::string directoryDisk;
	for (int i=0; i<rows; i++)
	{
		directorySS << directory << "/constant/polyMesh/sets-init/disk" << i+1 << "plusX";
		directoryDisk = directorySS.str();
		directorySS.str("");
		averageVelocity(dummyVel[i*4+0],directoryDisk,directoryU,directoryVol);
		directorySS << directory << "/constant/polyMesh/sets-init/disk" << i+1 << "minusX";
		directoryDisk = directorySS.str();
		directorySS.str("");
		averageVelocity(dummyVel[i*4+1],directoryDisk,directoryU,directoryVol);
		directorySS << directory << "/constant/polyMesh/sets-init/disk" << i+1 << "plusY";
		directoryDisk = directorySS.str();
		directorySS.str("");
		averageVelocity(dummyVel[i*4+2],directoryDisk,directoryU,directoryVol);
		directorySS << directory << "/constant/polyMesh/sets-init/disk" << i+1 << "minusY";
		directoryDisk = directorySS.str();
		directorySS.str("");
		averageVelocity(dummyVel[i*4+3],directoryDisk,directoryU,directoryVol);
	}
}

void ReadResults::dummyRotorVolumes(std::vector<double>& dummyVol, int finalIter, std::string directory) {
	int rows =  dummyVol.size()/4;
	std::ostringstream directorySS;
	directorySS << directory << "/" << finalIter << "/Vol";
	std::string directoryVol = directorySS.str();
	directorySS.str("");
	std::string directoryDisk;
	for (int i=0; i<rows; i++)
	{
		directorySS << directory << "/constant/polyMesh/sets-init/disk" << i+1 << "plusX";
		directoryDisk = directorySS.str();
		directorySS.str("");
		totalVolume(dummyVol[i*4+0],directoryDisk,directoryVol);
		directorySS << directory << "/constant/polyMesh/sets-init/disk" << i+1 << "minusX";
		directoryDisk = directorySS.str();
		directorySS.str("");
		totalVolume(dummyVol[i*4+1],directoryDisk,directoryVol);
		directorySS << directory << "/constant/polyMesh/sets-init/disk" << i+1 << "plusY";
		directoryDisk = directorySS.str();
		directorySS.str("");
		totalVolume(dummyVol[i*4+2],directoryDisk,directoryVol);
		directorySS << directory << "/constant/polyMesh/sets-init/disk" << i+1 << "minusY";
		directoryDisk = directorySS.str();
		directorySS.str("");
		totalVolume(dummyVol[i*4+3],directoryDisk,directoryVol);
	}
}

void ReadResults::adjointVelocities(std::vector<double>& adjointVel, int finalIter, std::string directory) {
	int rows =  adjointVel.size()/8;
	std::ostringstream directorySS;
	directorySS << directory << "/" << finalIter << "/U_a";
	std::string directoryU = directorySS.str();
	directorySS.str("");
	directorySS << directory << "/" << finalIter << "/Vol";
	std::string directoryVol = directorySS.str();
	directorySS.str("");
	std::string directoryDisk;
	for (int i=0; i<rows; i++)
	{
		directorySS << directory << "/constant/polyMesh/sets/disk" << i+1 << "plusX";
		directoryDisk = directorySS.str();
		directorySS.str("");
		averageVelocity(adjointVel[i*8+0],directoryDisk,directoryU,directoryVol);
		directorySS << directory << "/constant/polyMesh/sets/disk" << i+1 << "plusXX";
		directoryDisk = directorySS.str();
		directorySS.str("");
		averageVelocity(adjointVel[i*8+1],directoryDisk,directoryU,directoryVol);
		directorySS << directory << "/constant/polyMesh/sets/disk" << i+1 << "minusX";
		directoryDisk = directorySS.str();
		directorySS.str("");
		averageVelocity(adjointVel[i*8+2],directoryDisk,directoryU,directoryVol);
		directorySS << directory << "/constant/polyMesh/sets/disk" << i+1 << "minusXX";
		directoryDisk = directorySS.str();
		directorySS.str("");
		averageVelocity(adjointVel[i*8+3],directoryDisk,directoryU,directoryVol);
		directorySS << directory << "/constant/polyMesh/sets/disk" << i+1 << "plusY";
		directoryDisk = directorySS.str();
		directorySS.str("");
		averageVelocity(adjointVel[i*8+4],directoryDisk,directoryU,directoryVol);
		directorySS << directory << "/constant/polyMesh/sets/disk" << i+1 << "plusYY";
		directoryDisk = directorySS.str();
		directorySS.str("");
		averageVelocity(adjointVel[i*8+5],directoryDisk,directoryU,directoryVol);
		directorySS << directory << "/constant/polyMesh/sets/disk" << i+1 << "minusY";
		directoryDisk = directorySS.str();
		directorySS.str("");
		averageVelocity(adjointVel[i*8+6],directoryDisk,directoryU,directoryVol);
		directorySS << directory << "/constant/polyMesh/sets/disk" << i+1 << "minusYY";
		directoryDisk = directorySS.str();
		directorySS.str("");
		averageVelocity(adjointVel[i*8+7],directoryDisk,directoryU,directoryVol);
	}
}

void ReadResults::adjointVolumes(std::vector<double>& adjointVol, int finalIter, std::string directory) {
	int rows =  adjointVol.size()/8;
	std::ostringstream directorySS;
	directorySS << directory << "/" << finalIter << "/Vol";
	std::string directoryVol = directorySS.str();
	directorySS.str("");
	std::string directoryDisk;
	for (int i=0; i<rows; i++)
	{
		directorySS << directory << "/constant/polyMesh/sets/disk" << i+1 << "plusX";
		directoryDisk = directorySS.str();
		directorySS.str("");
		totalVolume(adjointVol[i*8+0],directoryDisk,directoryVol);
		directorySS << directory << "/constant/polyMesh/sets/disk" << i+1 << "plusXX";
		directoryDisk = directorySS.str();
		directorySS.str("");
		totalVolume(adjointVol[i*8+1],directoryDisk,directoryVol);
		directorySS << directory << "/constant/polyMesh/sets/disk" << i+1 << "minusX";
		directoryDisk = directorySS.str();
		directorySS.str("");
		totalVolume(adjointVol[i*8+2],directoryDisk,directoryVol);
		directorySS << directory << "/constant/polyMesh/sets/disk" << i+1 << "minusXX";
		directoryDisk = directorySS.str();
		directorySS.str("");
		totalVolume(adjointVol[i*8+3],directoryDisk,directoryVol);
		directorySS << directory << "/constant/polyMesh/sets/disk" << i+1 << "plusY";
		directoryDisk = directorySS.str();
		directorySS.str("");
		totalVolume(adjointVol[i*8+4],directoryDisk,directoryVol);
		directorySS << directory << "/constant/polyMesh/sets/disk" << i+1 << "plusYY";
		directoryDisk = directorySS.str();
		directorySS.str("");
		totalVolume(adjointVol[i*8+5],directoryDisk,directoryVol);
		directorySS << directory << "/constant/polyMesh/sets/disk" << i+1 << "minusY";
		directoryDisk = directorySS.str();
		directorySS.str("");
		totalVolume(adjointVol[i*8+6],directoryDisk,directoryVol);
		directorySS << directory << "/constant/polyMesh/sets/disk" << i+1 << "minusYY";
		directoryDisk = directorySS.str();
		directorySS.str("");
		totalVolume(adjointVol[i*8+7],directoryDisk,directoryVol);
	}
}

void ReadResults::averageVelocity(double& averageVel, std::string directoryDisk, std::string directoryU, std::string directoryVol) {
	std::ifstream ifs (directoryDisk.c_str());
	std::string line;
	for (int j=1; j<=19; j++) {
		std::getline(ifs, line);
	}
	std::istringstream tempSS(line);
	double totCells;
	tempSS >> totCells;
	tempSS.clear();
	std::vector<double> cell(totCells);
	std::getline(ifs, line);
	for (int j=0; j<totCells; j++) {
		std::getline(ifs, line);
		tempSS.str(line);
		tempSS >> cell[j];
		tempSS.clear();
	}
	ifs.close();
	std::sort(cell.begin(),cell.end());

	std::vector<double> vel(totCells);
	ifs.open(directoryU.c_str());
	for (int k=1; k<=22; k++) {
		std::getline(ifs, line);
	}
	int k=0;
	int j=0;
	while (std::getline(ifs, line) && j<totCells) {
		if (k == cell[j]) {
			tempSS.str(line);
			tempSS.seekg(2);
			tempSS >> vel[j];
			tempSS >> vel[j];
			tempSS.clear();
			j = j+1;
		}
		k = k+1;
	}
	ifs.close();

	std::vector<double> vol(totCells);
	ifs.open(directoryVol.c_str());
	for (int k=1; k<=20; k++) {
		std::getline(ifs, line);
	}
	k=0;
	j=0;
	line.erase(0,16);
	std::string lineTemp = "u";
	if (line[0] == lineTemp[0]) { //uniform
		line.erase(0,8);
		while (j<totCells) {
			tempSS.str(line);
			tempSS >> vol[j];
			tempSS.clear();
			j = j+1;
		}
	}
	else {
		std::getline(ifs, line);
		std::getline(ifs, line);
		while (std::getline(ifs, line) && j<totCells) {
			if (k == cell[j]) {
				tempSS.str(line);
				tempSS >> vol[j];
				tempSS.clear();
				j = j+1;
			}
			k = k+1;
		}
	}
	ifs.close();

	averageVel = 0;
	double totVol = 0;
	for (int j=0; j<totCells; j++) {
		averageVel = averageVel - vol[j]*vel[j];
		totVol = totVol + vol[j];
	}
	averageVel = averageVel/totVol;
}

void ReadResults::totalVolume(double& totVol, std::string directoryDisk, std::string directoryVol) {
	std::ifstream ifs (directoryDisk.c_str());
	std::string line;
	for (int j=1; j<=19; j++) {
		std::getline(ifs, line);
	}
	std::istringstream tempSS(line);
	double totCells;
	tempSS >> totCells;
	tempSS.clear();
	std::vector<double> cell(totCells);
	std::getline(ifs, line);
	for (int j=0; j<totCells; j++) {
		std::getline(ifs, line);
		tempSS.str(line);
		tempSS >> cell[j];
		tempSS.clear();
	}
	ifs.close();
	std::sort(cell.begin(),cell.end());

	std::vector<double> vol(totCells);
	ifs.open(directoryVol.c_str());
	for (int k=1; k<=20; k++) {
		std::getline(ifs, line);
	}
	int k=0;
	int j=0;
	line.erase(0,16);
	std::string lineTemp = "u";
	if (line[0] == lineTemp[0]) { //uniform
		line.erase(0,8);
		while (j<totCells) {
			tempSS.str(line);
			tempSS >> vol[j];
			tempSS.clear();
			j = j+1;
		}
	}
	else {
		std::getline(ifs, line);
		std::getline(ifs, line);
		while (std::getline(ifs, line) && j<totCells) {
			if (k == cell[j]) {
				tempSS.str(line);
				tempSS >> vol[j];
				tempSS.clear();
				j = j+1;
			}
			k = k+1;
		}
	}
	ifs.close();

	totVol = 0;
	for (int j=0; j<totCells; j++) {
		totVol = totVol + vol[j];
	}
}
