#include <math.h>
#include <cmath>
#include <vector>
#include <sys/stat.h>
#include <sys/types.h>
#include <sstream>
#include <iostream>
#include <fstream>
#include <string>

#include "ReadInputs.h"
#include "Geometry.h"
#include "DiffusionProperties.h"
#include "Simulation.h"
#include "BCvelocity.h"
#include "BCpressure.h"
#include "BCk.h"
#include "BCepsilon.h"
#include "BComega.h"
#include "BCnut.h"
#include "ReadResults.h"

#define pi 3.1416

int main() {

	ReadInputs inputs;
	Geometry mesh;
	DiffusionProperties diffusion;
	Simulation simulation;
	BCvelocity velocity;
	BCpressure pressure;
	BCk k;
	BCepsilon epsilon;
	BComega omega;
	BCnut nut;
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
	std::vector<double> Ustar(wsbinmax); // Friction velocity
	std::vector<double> Utop(wsbinmax); // Top velocity
	std::vector<double> ABLk(wsbinmax); // Turbulence kinetic energy
	std::vector<double> ABLepsilon(wsbinmax); // Turbulence dissipation rate
	std::vector<double> ABLepsilonTop(wsbinmax); // Top turbulence dissipation rate
	std::vector<double> ABLomega(wsbinmax); // Specific dissipation rate
	std::vector<double> ABLomegaTop(wsbinmax); // Top specific dissipation rate
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
		Ustar[i] = windSpeed[i]*kappa/log((h+h0+z0)/z0);
		Utop[i] = Ustar[i]/kappa*log((lz+z0)/z0);
		ABLk[i] = pow(Ustar[i],2)/pow(Cmu,0.5);
		ABLepsilon[i] = pow(Ustar[i],3)/(kappa*(h+h0+z0));
		ABLepsilonTop[i] = pow(Ustar[i],3)/(kappa*(lz+z0));
		ABLomega[i] = Ustar[i]/(kappa*pow(Cmu,0.5)*(h+h0+z0));
		ABLomegaTop[i] = Ustar[i]/(kappa*pow(Cmu,0.5)*(lz+z0));
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

	std::vector< std::vector< std::vector<double> > > rotorVel(wsbinmax,std::vector< std::vector<double> >(nwd,std::vector<double>(nt)));

	std::vector<double> Gradient(nt*2);

	std::ofstream resultFile;

	for (int i=0; i<nwd; i++)
	{
		std::vector< std::vector<double> > relPositions(nt,std::vector<double>(2));
		std::vector< std::vector<double> > tempRelPositions(nt,std::vector<double>(2));

		for (int z=0; z<nt; z++)
		{
			relPositions[z][0] = positions[z][0]*cos(windDir[i]*pi/180)-positions[z][1]*sin(windDir[i]*pi/180);
			relPositions[z][1] = positions[z][0]*sin(windDir[i]*pi/180)+positions[z][1]*cos(windDir[i]*pi/180);
		}
		for (int j=0; j<wsbinmax; j++)
		{
			if (windRose[j][i] > 1e-5)
			{
				resultFile.open("Results-CDGradient.txt",std::ofstream::app);
				resultFile << "Wind direction " << windDir[i] << " deg, wind speed " << windSpeed[j] << "\n\n";
				resultFile.close();
				std::vector<double> tempGradient(nt*2);
				std::vector<double> absGradient(nt*2);
				for (int z=0; z<nt*4; z++)
				{
					tempRelPositions = relPositions;
					std::ostringstream directorySS;
					resultFile.open("Results-CDGradient.txt",std::ofstream::app);
					if (z==0 || z%4 == 0)
					{
						directorySS << "./wd-" << windDir[i] << "-ws-" << windSpeed[j]-0.5 << "_5-Turb" << floor(z/4)+1 << "plusX";
						tempRelPositions[floor(z/4)][0] = tempRelPositions[floor(z/4)][0] + variation;
						resultFile << "Turbine " << floor(z/4)+1 << " plus delta x\n";
					}
					else if (z%4 == 1)
					{
						directorySS << "./wd-" << windDir[i] << "-ws-" << windSpeed[j]-0.5 << "_5-Turb" << floor(z/4)+1 << "minusX";
						tempRelPositions[floor(z/4)][0] = tempRelPositions[floor(z/4)][0] - variation;
						resultFile << "Turbine " << floor(z/4)+1 << " minus delta x\n";
					}
					else if (z%4 == 2)
					{
						directorySS << "./wd-" << windDir[i] << "-ws-" << windSpeed[j]-0.5 << "_5-Turb" << floor(z/4)+1 << "plusY";
						tempRelPositions[floor(z/4)][1] = tempRelPositions[floor(z/4)][1] + variation;
						resultFile << "Turbine " << floor(z/4)+1 << " plus delta y\n";
					}
					else if (z%4 == 3)
					{
						directorySS << "./wd-" << windDir[i] << "-ws-" << windSpeed[j]-0.5 << "_5-Turb" << floor(z/4)+1 << "minusY";
						tempRelPositions[floor(z/4)][1] = tempRelPositions[floor(z/4)][1] - variation;
						resultFile << "Turbine " << floor(z/4)+1 << " minus delta y\n";
					}
					resultFile.close();

					std::string directoryS = directorySS.str();
					mkdir(directoryS.c_str(), S_IRWXU | S_IRWXG | S_IROTH | S_IXOTH);
					directoryS = directorySS.str()+"/0";
					mkdir(directoryS.c_str(), S_IRWXU | S_IRWXG | S_IROTH | S_IXOTH);
					directoryS = directorySS.str()+"/0.primal";
					mkdir(directoryS.c_str(), S_IRWXU | S_IRWXG | S_IROTH | S_IXOTH);
					directoryS = directorySS.str()+"/constant";
					mkdir(directoryS.c_str(), S_IRWXU | S_IRWXG | S_IROTH | S_IXOTH);
					directoryS = directorySS.str()+"/constant/polyMesh";
					mkdir(directoryS.c_str(), S_IRWXU | S_IRWXG | S_IROTH | S_IXOTH);
					directoryS = directorySS.str()+"/system";
					mkdir(directoryS.c_str(), S_IRWXU | S_IRWXG | S_IROTH | S_IXOTH);

					if (dim == 2) 
					{
						mesh.twodLayout(lx,ly,delta,directorySS.str());
						pressure.twodLayout(directorySS.str());
						velocity.twodLayout(windSpeed[j],directorySS.str());
						if (turb > 0) 
						{
							diffusion.turbulentKO(directorySS.str());
							k.twodLayout(ABLk[j],directorySS.str());
							omega.twodLayout(ABLomega[j],directorySS.str());
							epsilon.twodLayout(ABLepsilon[j],directorySS.str());
							nut.twodLayout(directorySS.str());
						}
					}
					else if (dim == 3)
					{
						if (terrain == 0) 
						{
							mesh.threedLayout(lx,ly,lz,d,delta,directorySS.str());
							mesh.threedLayoutRefine(tempRelPositions,d,h,directorySS.str());
						}
						else
						{
							mesh.remapTerrain(windDir[i],directorySS.str());
							mesh.remapTurbines(tempRelPositions,directorySS.str());
							mesh.threedLayoutComplexTerrain(lx,ly,lz,d,delta,directorySS.str());
							mesh.threedLayoutRefineComplexTerrain(d,h,directorySS.str());
						}
						pressure.threedLayout(directorySS.str());
						velocity.threedLayout(windSpeed[j],Utop[j],lz,z0,directorySS.str());
						diffusion.turbulentKO(directorySS.str());
						k.threedLayout(ABLk[j],directorySS.str());
						omega.threedLayout(ABLomega[j],ABLomegaTop[j],Ustar[j],z0,directorySS.str());
						epsilon.threedLayout(ABLepsilon[j],ABLepsilonTop[j],Utop[j],lz,z0,directorySS.str());
						nut.threedLayout(z0,directorySS.str());
					}
					diffusion.viscosity(nu,directorySS.str());
					simulation.schemes(directorySS.str());
					simulation.solvers(directorySS.str());
					simulation.options(nt,T[j],directorySS.str());
					if (terrain == 0) 
					{
						simulation.topoSet(tempRelPositions,d,h,variation,dim,directorySS.str());
					}
					else
					{
						simulation.topoSetComplexTerrain(d,h,variation,directorySS.str());
					}
					simulation.control(nIter,turb,directorySS.str());
					simulation.run(np,turb,dim,directorySS.str());
					results.finalIter(finalIter,turb,directorySS.str());
					results.rotorVelocities(rotorVel[j][i],finalIter,directorySS.str());
					for (int k=0; k<nt; k++)
					{
						if (z == 0 || z%4 == 0)
						{
							tempGradient[floor(z/4)*2+0] = tempGradient[floor(z/4)*2+0] + 8760*T[j]*rotorVel[j][i][k]/(2*variation);
						}
						else if (z%4 == 1)
						{
							tempGradient[floor(z/4)*2+0] = tempGradient[floor(z/4)*2+0] - 8760*T[j]*rotorVel[j][i][k]/(2*variation);
						}
						else if (z%4 == 2)
						{
							tempGradient[floor(z/4)*2+1] = tempGradient[floor(z/4)*2+1] + 8760*T[j]*rotorVel[j][i][k]/(2*variation);
						}
						else if (z%4 == 3)
						{
							tempGradient[floor(z/4)*2+1] = tempGradient[floor(z/4)*2+1] - 8760*T[j]*rotorVel[j][i][k]/(2*variation);
						}
					}
					simulation.cleanPrimal(directorySS.str());
					directorySS.str("");
					resultFile.open("Results-CDGradient.txt",std::ofstream::app);
					resultFile << "Finished after " << finalIter << " iterations\n";
					resultFile << "Thrust " << T[j] << " N\n";
					for (int k=0; k<nt; k++)
					{
						resultFile << "Wind turbine " << k+1 << ", average rotor speed " << rotorVel[j][i][k] << " m/s\n";
					}
					resultFile << "\n";
					resultFile.close();
				}
				for (int k=0; k<nt; k++)
				{
					absGradient[k*2+0] = tempGradient[k*2+0]*cos(-windDir[i]*pi/180)-tempGradient[k*2+1]*sin(-windDir[i]*pi/180);
					absGradient[k*2+1] = tempGradient[k*2+0]*sin(-windDir[i]*pi/180)+tempGradient[k*2+1]*cos(-windDir[i]*pi/180);
					Gradient[k*2+0] = Gradient[k*2+0] + windRose[j][i]*absGradient[k*2+0];
					Gradient[k*2+1] = Gradient[k*2+1] + windRose[j][i]*absGradient[k*2+1];
				}
				resultFile.open("Results-CDGradient.txt",std::ofstream::app);
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
	resultFile.open("Results-CDGradient.txt",std::ofstream::app);
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
	resultFile << "\n";
	resultFile.close();

	return 0;
}
