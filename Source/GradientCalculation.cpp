#include <math.h>
#include <cmath>
#include <vector>
#include <sys/stat.h>
#include <sys/types.h>
#include <sstream>
#include <iostream>
#include <fstream>
#include <string>
#include <algorithm>

#include "ReadInputs.h"
#include "DiffusionProperties.h"
#include "Simulation.h"
#include "BCvelocity.h"
#include "BCpressure.h"
#include "BCk.h"
#include "BCepsilon.h"
#include "BComega.h"
#include "ReadResults.h"

#define pi 3.1416

int main() {

	ReadInputs inputs;
	DiffusionProperties diffusion;
	Simulation simulation;
	BCvelocity velocity;
	BCpressure pressure;
	BCk k;
	BCepsilon epsilon;
	BComega omega;
	ReadResults results;

	int nt; // Number of wind turbines
	int nwd; // Number of wind directions
	double d; // Turbine diameter
	double h; // Hub height
	double h0; // Turbine base height
	double z0; // Surface roughness
	double lx; // Streamwise lenght
	double ly; // Lateral length
	double lz; // Domain height
	double delta; // Mesh resolution
	int turb; // Diffusion type
	int dim; // Problem dimensions
	int np; // Number of processor
	int terrain; // Terrain type
	double height; // Height of Gaussian hill
	double spread; // Spread of Gaussian hill

	inputs.setup(d,h,h0,z0,lx,ly,lz,delta,turb,dim,np,terrain,height,spread);

	int wsbinmax = 30; // Maximum wind speed bin

	std::vector<double> Ct(wsbinmax); // Thrust coefficient
	std::vector< std::vector<double> > positions; // Turbine positions
	std::vector<double> windDir; // Wind directions
	std::vector< std::vector<double> > windRose; // Wind rose frequency

	inputs.thrustCoeff(Ct);
	inputs.positions(nt,positions);
	inputs.windRose(wsbinmax,nwd,windDir,windRose);

	double totFreq = 0;
	for (int i=0; i<nwd; i++)
	{
		for (int j=0; j<wsbinmax; j++)
		{
			totFreq = totFreq+windRose[j][i];
		}
	}
	if (std::abs(totFreq-1)>0.01)
	{
		std::cout << "Warning: Total frequency is off from 1 by "<<totFreq-1<<"\n";
	}

	double rho = 1.2; // Air density
	double nu = 1.5e-05; // Air kinetic viscosity
	double A = pi*pow(d,2)/4.0; // Rotor spanned area
	
	double kappa = 0.41; // Von-Karman constant
	double Cmu = 0.0333; // Turbulent constant

	double variation; // Variation of the design variables
	if (dim == 2) 
	{
		variation = delta;
	}
	else if (dim == 3)
	{
		variation = delta/4;
	}

	std::vector<double> windSpeed(wsbinmax);
	std::vector<double> T(wsbinmax); // Thrust

	for (int i=0; i<wsbinmax; i++)
	{
		windSpeed[i] = i+0.5;
		if (dim == 2) 
		{
			T[i] = 0.5*rho*d*pow(windSpeed[i],2)*Ct[i];
		}
		else if (dim == 3)
		{
			T[i] = 0.5*rho*A*pow(windSpeed[i],2)*Ct[i];
		}
	}


	double dt = 1;
	int nIter = 10000;
	int finalIter = 0;

	if (turb == 0) 
	{
		for (int i=0; i<wsbinmax; i++)
		{
			T[i] = 0.1*T[i];
		}
	}

	std::vector< std::vector< std::vector<double> > > adjointVel(wsbinmax,std::vector< std::vector<double> >(nwd,std::vector<double>(nt*8)));
	std::vector< std::vector< std::vector<double> > > adjointVol(wsbinmax,std::vector< std::vector<double> >(nwd,std::vector<double>(nt*8)));
	std::vector< std::vector< std::vector<double> > > dummyVel(wsbinmax,std::vector< std::vector<double> >(nwd,std::vector<double>(nt*4)));
	std::vector< std::vector< std::vector<double> > > dummyVol(wsbinmax,std::vector< std::vector<double> >(nwd,std::vector<double>(nt*4)));
	std::vector< std::vector< std::vector<double> > > rotorVol(wsbinmax,std::vector< std::vector<double> >(nwd,std::vector<double>(nt)));

	std::vector< std::vector< std::vector<double> > > rotorVel(wsbinmax,std::vector< std::vector<double> >(nwd,std::vector<double>(nt)));

	std::vector<double> Gradient(nt*2);

	std::ofstream resultFile;

	for (int i=0; i<nwd; i++)
	{
		std::vector< std::vector<double> > relPositions(nt,std::vector<double>(2));
		for (int z=0; z<nt; z++)
		{
			relPositions[z][0] = positions[z][0]*cos(windDir[i]*pi/180)-positions[z][1]*sin(windDir[i]*pi/180);
			relPositions[z][1] = positions[z][0]*sin(windDir[i]*pi/180)+positions[z][1]*cos(windDir[i]*pi/180);
		}
		for (int j=0; j<wsbinmax; j++)
		{
			if (windRose[j][i] > 1e-5)
			{
				double check = 0;
				std::vector<double> tempGradient(nt*2);
				std::ostringstream directorySS;
				while (check < 1 || isnan(check) || isinf(check) || finalIter == 0)
				{
					directorySS.str("");
					directorySS << "./wd-" << windDir[i] << "-ws-" << windSpeed[j]-0.5 << "_5";
					std::string command = "cd " + directorySS.str() + "; ls -d *Output | xargs -I {} rm {}";
					system(command.c_str());
					command = "cd " + directorySS.str() + "; ls -d 1* | xargs -I {} rm -r {}";
					system(command.c_str());
					command = "cd " + directorySS.str() + "; ls -d 2* | xargs -I {} rm -r {}";
					system(command.c_str());
					command = "cd " + directorySS.str() + "; ls -d 3* | xargs -I {} rm -r {}";
					system(command.c_str());
					command = "cd " + directorySS.str() + "; ls -d 4* | xargs -I {} rm -r {}";
					system(command.c_str());
					command = "cd " + directorySS.str() + "; ls -d 5* | xargs -I {} rm -r {}";
					system(command.c_str());
					command = "cd " + directorySS.str() + "; ls -d 6* | xargs -I {} rm -r {}";
					system(command.c_str());
					command = "cd " + directorySS.str() + "; ls -d 7* | xargs -I {} rm -r {}";
					system(command.c_str());
					command = "cd " + directorySS.str() + "; ls -d 8* | xargs -I {} rm -r {}";
					system(command.c_str());
					command = "cd " + directorySS.str() + "; ls -d 9* | xargs -I {} rm -r {}";
					system(command.c_str());
					std::string directoryS = directorySS.str()+"/0";
					command = "rm -r " + directoryS;
					system(command.c_str());
					mkdir(directoryS.c_str(), S_IRWXU | S_IRWXG | S_IROTH | S_IXOTH);
					directoryS = directorySS.str()+"/constant";
					command = "rm -r " + directoryS;
					system(command.c_str());
					mkdir(directoryS.c_str(), S_IRWXU | S_IRWXG | S_IROTH | S_IXOTH);
					directoryS = directorySS.str()+"/system";
					command = "rm -r " + directoryS;
					system(command.c_str());
					mkdir(directoryS.c_str(), S_IRWXU | S_IRWXG | S_IROTH | S_IXOTH);
	
					if (dim == 2) 
					{
						pressure.a_twodLayout(directorySS.str());
						velocity.a_twodLayout(directorySS.str());
						k.a_twodLayout(directorySS.str());
						epsilon.a_twodLayout(directorySS.str());
						omega.a_twodLayout(directorySS.str());
					}
					else if (dim == 3)
					{
						pressure.a_threedLayout(directorySS.str());
						velocity.a_threedLayout(directorySS.str());
						k.a_threedLayout(directorySS.str());
						epsilon.a_threedLayout(directorySS.str());
						omega.a_threedLayout(directorySS.str());
					}
					simulation.a_schemes(directorySS.str());
					simulation.a_solvers(turb,directorySS.str());
					simulation.a_options(nt,T[j],directorySS.str());
					if (terrain == 0) 
					{
						simulation.topoSet(relPositions,d,h,variation,dim,directorySS.str());
						simulation.a_topoSet(relPositions,lx,ly,lz,d,h,variation,dim,directorySS.str());
					}
					else
					{
						simulation.topoSetComplexTerrain(d,h,variation,directorySS.str());
						simulation.a_topoSetCompleTerrain(lx,ly,lz,d,h,variation,dim,directorySS.str());
					}
					simulation.a_control(nIter,turb,directorySS.str());
					simulation.a_run(np,turb,dim,directorySS.str());
					results.a_finalIter(finalIter,turb,directorySS.str());
					results.adjointVelocities(adjointVel[j][i],finalIter,directorySS.str());
					results.adjointVolumes(adjointVol[j][i],finalIter,directorySS.str());
					results.dummyRotorVelocities(dummyVel[j][i],finalIter,directorySS.str());
					results.dummyRotorVolumes(dummyVol[j][i],finalIter,directorySS.str());
					results.rotorVolumes(rotorVol[j][i],finalIter,directorySS.str());
					results.rotorVelocities(rotorVel[j][i],finalIter,directorySS.str());
					check = 0;
					for (int k=0; k<nt*2; k++)
					{
						tempGradient[k] = 0;
					}
					for (int k=0; k<nt; k++)
					{
						tempGradient[k*2+0] = tempGradient[k*2+0] + 8760*T[j]/rotorVol[j][i][k]*dummyVel[j][i][k*4+0]*dummyVol[j][i][k*4+0]/(2*variation);
						tempGradient[k*2+0] = tempGradient[k*2+0] - 8760*T[j]/rotorVol[j][i][k]*dummyVel[j][i][k*4+1]*dummyVol[j][i][k*4+1]/(2*variation);
						tempGradient[k*2+1] = tempGradient[k*2+1] + 8760*T[j]/rotorVol[j][i][k]*dummyVel[j][i][k*4+2]*dummyVol[j][i][k*4+2]/(2*variation);
						tempGradient[k*2+1] = tempGradient[k*2+1] - 8760*T[j]/rotorVol[j][i][k]*dummyVel[j][i][k*4+3]*dummyVol[j][i][k*4+3]/(2*variation);
						tempGradient[k*2+0] = tempGradient[k*2+0] + 8760*T[j]/rotorVol[j][i][k]*adjointVel[j][i][k*8+0]*adjointVol[j][i][k*8+0]/(2*variation);
						tempGradient[k*2+0] = tempGradient[k*2+0] - 8760*T[j]/rotorVol[j][i][k]*adjointVel[j][i][k*8+1]*adjointVol[j][i][k*8+1]/(2*variation);
						tempGradient[k*2+0] = tempGradient[k*2+0] + 8760*T[j]/rotorVol[j][i][k]*adjointVel[j][i][k*8+2]*adjointVol[j][i][k*8+2]/(2*variation);
						tempGradient[k*2+0] = tempGradient[k*2+0] - 8760*T[j]/rotorVol[j][i][k]*adjointVel[j][i][k*8+3]*adjointVol[j][i][k*8+3]/(2*variation);
						tempGradient[k*2+1] = tempGradient[k*2+1] + 8760*T[j]/rotorVol[j][i][k]*adjointVel[j][i][k*8+4]*adjointVol[j][i][k*8+4]/(2*variation);
						tempGradient[k*2+1] = tempGradient[k*2+1] - 8760*T[j]/rotorVol[j][i][k]*adjointVel[j][i][k*8+5]*adjointVol[j][i][k*8+5]/(2*variation);
						tempGradient[k*2+1] = tempGradient[k*2+1] + 8760*T[j]/rotorVol[j][i][k]*adjointVel[j][i][k*8+6]*adjointVol[j][i][k*8+6]/(2*variation);
						tempGradient[k*2+1] = tempGradient[k*2+1] - 8760*T[j]/rotorVol[j][i][k]*adjointVel[j][i][k*8+7]*adjointVol[j][i][k*8+7]/(2*variation);
						check = check + abs(tempGradient[k*2+0]) + abs(tempGradient[k*2+1]);
					}
				}
				std::vector<double> absGradient(nt*2);
				for (int k=0; k<nt; k++)
				{
					//if (rotorVel[j][i][k] > 0.95*(*std::max_element(rotorVel[j][i].begin(), rotorVel[j][i].end()))) {
					//	tempGradient[k*2+1] = 0;
					//}
					absGradient[k*2+0] = tempGradient[k*2+0]*cos(-windDir[i]*pi/180)-tempGradient[k*2+1]*sin(-windDir[i]*pi/180);
					absGradient[k*2+1] = tempGradient[k*2+0]*sin(-windDir[i]*pi/180)+tempGradient[k*2+1]*cos(-windDir[i]*pi/180);
					Gradient[k*2+0] = Gradient[k*2+0] + windRose[j][i]*absGradient[k*2+0];
					Gradient[k*2+1] = Gradient[k*2+1] + windRose[j][i]*absGradient[k*2+1];
				}
				simulation.cleanAdjoint(directorySS.str());
				directorySS.str("");
				resultFile.open("Results-Gradient.txt",std::ofstream::app);
				resultFile << "Wind direction " << windDir[i] << " deg, wind speed " << windSpeed[j] << " m/s\n\n";
				resultFile << "Finished after " << finalIter << " iterations\n";
				resultFile << "The gradient is\n";
				for (int k=0; k<nt; k++)
				{
					resultFile << tempGradient[k*2+0]/1000 << " " << tempGradient[k*2+1]/1000 << "\n";
				}
				resultFile << "kWh/m\nor\n";
				for (int k=0; k<nt; k++)
				{
					resultFile << sqrt(pow(tempGradient[k*2+0]/1000,2) + pow(tempGradient[k*2+1]/1000,2)) << " kWh/m, " << atan2(tempGradient[k*2+1]/1000,tempGradient[k*2+0]/1000)*180/pi << " deg\n";
				}
				resultFile << "\n";
				resultFile.close();
			}
		}
	}
	resultFile.open("Results-Gradient.txt",std::ofstream::app);
	resultFile << "The total gradient is\n";
	for (int k=0; k<nt; k++)
	{
		resultFile << Gradient[k*2+0]/1000 << " " << Gradient[k*2+1]/1000 << "\n";
	}
	resultFile << "kWh/m\nor\n";
	for (int k=0; k<nt; k++)
	{
		resultFile << sqrt(pow(Gradient[k*2+0]/1000,2) + pow(Gradient[k*2+1]/1000,2)) << " kWh/m, " << atan2(Gradient[k*2+1]/1000,Gradient[k*2+0]/1000)*180/pi << " deg\n";
	}
	resultFile.close();

	return 0;
}
