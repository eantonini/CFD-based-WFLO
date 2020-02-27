#include "Geometry.h"
#include <fstream>
#include <sstream>
#include <math.h>
#include <vector>
#include <string>
#include <stdlib.h>
#include <stdio.h>
#define pi 3.1416

Geometry::Geometry() {
}

void Geometry::twodLayout(double lx, double ly, double delta, std::string directory) {
	std::string directoryS = directory+"/system/blockMeshDict";
	std::ofstream ofs;
	ofs.open (directoryS.c_str());
	ofs << "FoamFile{version 2.0; format ascii; class dictionary; location \"system\"; object blockMeshDict;}\n";
	ofs << "convertToMeters 1; vertices (("<<-lx/2<<" "<<-ly/2<<" 0) ("<<lx/2<<" "<<-ly/2<<" 0) ("<<lx/2<<" "<<ly/2<<" 0) ("<<-lx/2<<" "<<ly/2<<" 0)\n";
	ofs << "("<<-lx/2<<" "<<-ly/2<<" 1) ("<<lx/2<<" "<<-ly/2<<" 1) ("<<lx/2<<" "<<ly/2<<" 1) ("<<-lx/2<<" "<<ly/2<<" 1));\n";
	ofs << "blocks (hex (0 1 2 3 4 5 6 7) ("<<ceil(lx/delta)<<" "<<ceil(ly/delta)<<" 1) simpleGrading (1 1 1));\n";
	ofs << "edges ();\n";
	ofs << "boundary (\n";
	ofs << "inlet{type patch; faces ((3 2 6 7));}\n";
	ofs << "outlet{type patch; faces ((0 1 5 4));}\n";
	ofs << "upAndDown{type empty; faces ((0 1 2 3) (4 5 6 7));}\n";
	ofs << "sideE{type patch; faces ((1 5 6 2));}\n";
	ofs << "sideW{type patch; faces ((0 4 7 3));});\n";
	ofs << "mergePatchPairs ();";
	ofs.close();
}

void Geometry::threedLayout(double lx, double ly, double lz, double d, double delta, std::string directory) {
	int nz;
	int kz;
	if (delta/d < 0.15 && delta/d >= 0.08) {
		nz = 22;
		kz = 64;
	}
	else if (delta/d < 0.08 && delta/d >= 0.04) {
		nz = 36;
		kz = 36;
	}
	std::string directoryS = directory+"/system/blockMeshDict";
	std::ofstream ofs;
	ofs.open (directoryS.c_str());
	ofs << "FoamFile{version 2.0; format ascii; class dictionary; location \"system\"; object blockMeshDict;}\n";
	ofs << "convertToMeters 1; vertices (("<<-lx/2<<" "<<-ly/2<<" 0) ("<<lx/2<<" "<<-ly/2<<" 0) ("<<lx/2<<" "<<ly/2<<" 0) ("<<-lx/2<<" "<<ly/2<<" 0)\n";
	ofs << "("<<-lx/2<<" "<<-ly/2<<" "<<d/2<<") ("<<lx/2<<" "<<-ly/2<<" "<<d/2<<") ("<<lx/2<<" "<<ly/2<<" "<<d/2<<") ("<<-lx/2<<" "<<ly/2<<" "<<d/2<<")\n";
	ofs << "("<<-lx/2<<" "<<-ly/2<<" "<<lz<<") ("<<lx/2<<" "<<-ly/2<<" "<<lz<<") ("<<lx/2<<" "<<ly/2<<" "<<lz<<") ("<<-lx/2<<" "<<ly/2<<" "<<lz<<"));\n";
	ofs << "blocks (hex (0 1 2 3 4 5 6 7) ("<<ceil(lx/delta)<<" "<<ceil(ly/delta)<<" "<<nz<<") simpleGrading (1 1 "<<kz<<")\n";
	ofs << "hex (4 5 6 7 8 9 10 11) ("<<ceil(lx/delta)<<" "<<ceil(ly/delta)<<" "<<ceil((lz-d/2)/delta)<<") simpleGrading (1 1 1));\n";
	ofs << "edges ();\n";
	ofs << "boundary (\n";
	ofs << "inlet{type patch; faces ((3 2 6 7) (7 6 10 11));}\n";
	ofs << "outlet{type patch; faces ((0 1 5 4) (4 5 9 8));}\n";
	ofs << "terrain{type wall; faces ((0 1 2 3));}\n";
	ofs << "top{type patch; faces ((8 9 10 11));}\n";
	ofs << "sideE{type patch; faces ((1 5 6 2) (5 9 10 6));}\n";
	ofs << "sideW{type patch; faces ((0 4 7 3) (4 8 11 7));});\n";
	ofs << "mergePatchPairs ();";
	ofs.close();
}

void Geometry::threedLayoutRefine(std::vector< std::vector<double> > positions, double d, double h, std::string directory) {
	int rows =  positions.size();
	std::string directoryS = directory+"/system/snappyHexMeshDict";
	std::ofstream ofs;
	ofs.open (directoryS.c_str());
	ofs << "FoamFile{version 2.0; format ascii; class dictionary; location \"system\"; object snappyHexMeshDict;}\n";
	ofs << "castellatedMesh true; snap false; addLayers false;\n";
	ofs << "geometry{\n";
	for (int i=0; i<rows; i++)
	{
		ofs << "innerCylinder"<<i+1<<"{type searchableCylinder; ";
		ofs << "point1 ("<<positions[i][0]<<" "<<positions[i][1]-d/5<<" "<<h<<"); point2 ("<<positions[i][0]<<" "<<positions[i][1]+d/5<<" "<<h<<"); radius "<<d/2*1.2<<";}\n";
		ofs << "outerCylinder"<<i+1<<"{type searchableCylinder; ";
		ofs << "point1 ("<<positions[i][0]<<" "<<positions[i][1]-d*1.5<<" "<<h<<"); point2 ("<<positions[i][0]<<" "<<positions[i][1]+d*1.5<<" "<<h<<"); radius "<<d/2*1.5<<";}\n";
	}
	ofs << "}\n";
	ofs << "castellatedMeshControls{maxLocalCells 10000000; maxGlobalCells 1000000000; minRefinementCells 0; maxLoadUnbalance 0.10; nCellsBetweenLevels 3; features ();\n";
	ofs << "refinementSurfaces{}\n";
	ofs << "resolveFeatureAngle 30; refinementRegions{\n";
	for (int i=0; i<rows; i++)
	{
		ofs << "innerCylinder"<<i+1<<"{mode inside; levels ((2 2));}\n";
		ofs << "outerCylinder"<<i+1<<"{mode inside; levels ((1 1));}\n";
	}
	ofs << "}\n";
	ofs << "locationInMesh (0.1 0.1 0.1); allowFreeStandingZoneFaces true;}\n";
	ofs << "snapControls{nSmoothPatch 4; tolerance 4; nSolveIter 30; nRelaxIter 10;}\n";
	ofs << "meshQualityControls{maxNonOrtho 65; maxBoundarySkewness 20; maxInternalSkewness 4; maxConcave 80; minVol 1e-13; minTetQuality 1e-15; minArea -1; minTwist 0.02;\n";
	ofs << "minDeterminant 0.001; minFaceWeight 0.05; minVolRatio 0.01; minTriangleTwist -1; nSmoothScale 4; errorReduction 0.75; relaxed{maxNonOrtho 75;}}\n";
	ofs << "addLayersControls{relativeSizes true; layers{} expansionRatio 1.0; finalLayerThickness 0.5; minThickness 0.25; nGrow 0; featureAngle 60; nRelaxIter 5;\n";
	ofs << "nSmoothSurfaceNormals 1; nSmoothNormals 3; nSmoothThickness 10; maxFaceThicknessRatio 0.5; maxThicknessToMedialRatio 0.3; minMedianAxisAngle 90;\n";
	ofs << "nBufferCellsNoExtrude 0; nLayerIter 50;}; mergeTolerance 1e-6;\n";
	ofs.close();
}

void Geometry::gaussianHill(double hight, double spread, double domainWidth, double resolution) {
	int n = domainWidth/resolution;
	double x0 = 0;
	double y0 = 0;
	double x;
	double y;
	double z;
	std::ofstream terrain ("terrain.txt");
	terrain << "X\tY\tZ";
	for (int i = 0; i < n+1; i++) {
		for (int j = 0; j < n+1; j++) {
			x = -domainWidth/2+i*resolution;
			y = -domainWidth/2+j*resolution;
			z = hight*exp(-(pow(x-x0,2)/(2*pow(spread,2))+pow(y-y0,2)/(2*pow(spread,2))));
			terrain << "\n" << x << "\t" << y << "\t" << z;
		}
	}
	terrain.close();
}

void Geometry::remapTerrain(double windDir, std::string directory) {
	std::vector< std::vector<double> > nodes;
	std::vector<double> temp(3);
	std::string line;
	std::istringstream iss;
	std::ifstream inputs ("terrain.txt");
	std::getline(inputs, line);
	int n = 0;
	while (std::getline(inputs, line))
	{
		iss.str(line);
		iss >> temp[0];
		iss >> temp[1];
		iss >> temp[2];
		nodes.push_back(temp);
		n = n+1;
		iss.clear();
	}
	inputs.close();
	double relX;
	double relY;
	double relZ;
	std::string directoryS = directory+"/terrainRotated.txt";
	std::ofstream terrain (directoryS.c_str());
	terrain << "X\tY\tZ";
	for (int i = 0; i < n; i++) {
		relX = nodes[i][0]*cos(windDir*pi/180)-nodes[i][1]*sin(windDir*pi/180);
		relY = nodes[i][0]*sin(windDir*pi/180)+nodes[i][1]*cos(windDir*pi/180);
		relZ = nodes[i][2];
		terrain << "\n" << relX << "\t" << relY << "\t" << relZ;
	}
	terrain.close();
	std::string str = "python ./RemapTerrain.py " + directory + "/";
	system(str.c_str());
}

void Geometry::remapTurbines(std::vector< std::vector<double> > positions, std::string directory) {
	int n =  positions.size();
	std::string directoryS = directory+"/turbinesRotated.txt";
	std::ofstream turbines (directoryS.c_str());
	turbines << "X\tY";
	for (int i = 0; i < n; i++) {
		turbines << "\n" << positions[i][0] << "\t" << positions[i][1];
	}
	turbines.close();
	std::string str = "python ./RemapTurbines.py " + directory + "/";
	system(str.c_str());
}

void Geometry::threedLayoutComplexTerrain(double lx, double ly, double lz, double d, double delta, std::string directory) {
	std::vector< std::vector<double> > nodes;
	std::vector<double> temp(3);
	std::string line;
	std::istringstream iss;
	std::string directoryS = directory+"/terrainMeshRotated.txt";
	std::ifstream inputs (directoryS.c_str());
	std::getline(inputs, line);
	int n = 0;
	while (std::getline(inputs, line))
	{
		iss.str(line);
		iss >> temp[0];
		iss >> temp[1];
		iss >> temp[2];
		nodes.push_back(temp);
		n = n+1;
		iss.clear();
	}
	inputs.close();
	int nz;
	int kz;
	if (delta/d < 0.15 && delta/d >= 0.08) {
		nz = 22;
		kz = 64;
	}
	else if (delta/d < 0.08 && delta/d >= 0.04) {
		nz = 36;
		kz = 36;
	}
	int nDiv = ceil(lx/delta);
	directoryS = directory+"/system/blockMeshDict";
	std::ofstream ofs;
	ofs.open (directoryS.c_str());
	ofs << "FoamFile{version 2.0; format ascii; class dictionary; location \"system\"; object blockMeshDict;}\n";
	ofs << "convertToMeters 1; vertices (\n";
	ofs << "("<<nodes[0][0]<<" "<<nodes[0][1]<<" "<<nodes[0][2]<<") ("<<nodes[1][0]<<" "<<nodes[1][1]<<" "<<nodes[1][2]<<") ("<<nodes[(nDiv+1)*nDiv+1][0]<<" "<<nodes[(nDiv+1)*nDiv+1][1]<<" "<<nodes[(nDiv+1)*nDiv+1][2]<<") ("<<nodes[(nDiv+1)*nDiv][0]<<" "<<nodes[(nDiv+1)*nDiv][1]<<" "<<nodes[(nDiv+1)*nDiv][2]<<")\n";
	ofs << "("<<nodes[0][0]<<" "<<nodes[0][1]<<" "<<nodes[0][2]+d/2<<") ("<<nodes[1][0]<<" "<<nodes[1][1]<<" "<<nodes[1][2]+d/2<<") ("<<nodes[(nDiv+1)*nDiv+1][0]<<" "<<nodes[(nDiv+1)*nDiv+1][1]<<" "<<nodes[(nDiv+1)*nDiv+1][2]+d/2<<") ("<<nodes[(nDiv+1)*nDiv][0]<<" "<<nodes[(nDiv+1)*nDiv][1]<<" "<<nodes[(nDiv+1)*nDiv][2]+d/2<<")\n";
	ofs << "("<<nodes[0][0]<<" "<<nodes[0][1]<<" "<<lz<<") ("<<nodes[1][0]<<" "<<nodes[1][1]<<" "<<lz<<") ("<<nodes[(nDiv+1)*nDiv+1][0]<<" "<<nodes[(nDiv+1)*nDiv+1][1]<<" "<<lz<<") ("<<nodes[(nDiv+1)*nDiv][0]<<" "<<nodes[(nDiv+1)*nDiv][1]<<" "<<lz<<")\n";
	for (int i = 2; i < nDiv+1; i++) {
		ofs << "("<<nodes[i][0]<<" "<<nodes[i][1]<<" "<<nodes[i][2]<<") ("<<nodes[(nDiv+1)*nDiv+i][0]<<" "<<nodes[(nDiv+1)*nDiv+i][1]<<" "<<nodes[(nDiv+1)*nDiv+i][2]<<") ";
		ofs << "("<<nodes[i][0]<<" "<<nodes[i][1]<<" "<<nodes[i][2]+d/2<<") ("<<nodes[(nDiv+1)*nDiv+i][0]<<" "<<nodes[(nDiv+1)*nDiv+i][1]<<" "<<nodes[(nDiv+1)*nDiv+i][2]+d/2<<") ";
		ofs << "("<<nodes[i][0]<<" "<<nodes[i][1]<<" "<<lz<<") ("<<nodes[(nDiv+1)*nDiv+i][0]<<" "<<nodes[(nDiv+1)*nDiv+i][1]<<" "<<lz<<")\n";
	}
	ofs << ");\n";
	ofs << "blocks (";
	ofs << "hex (0 1 2 3 4 5 6 7) (1 "<<nDiv<<" "<<nz<<") simpleGrading (1 1 "<<kz<<")\n";
	ofs << "hex (4 5 6 7 8 9 10 11) (1 "<<nDiv<<" "<<ceil((lz-d/2)/delta)<<") simpleGrading (1 1 1)\n";
	ofs << "hex (1 12 13 2 5 14 15 6) (1 "<<nDiv<<" "<<nz<<") simpleGrading (1 1 "<<kz<<")\n";
	ofs << "hex (5 14 15 6 9 16 17 10) (1 "<<nDiv<<" "<<ceil((lz-d/2)/delta)<<") simpleGrading (1 1 1)\n";
	for (int i = 2; i < nDiv; i++) {
		ofs << "hex ("<<i*6<<" "<<i*6+6<<" "<<i*6+7<<" "<<i*6+1<<" "<<i*6+2<<" "<<i*6+8<<" "<<i*6+9<<" "<<i*6+3<<") (1 "<<nDiv<<" "<<nz<<") simpleGrading (1 1 "<<kz<<")\n";
		ofs << "hex ("<<i*6+2<<" "<<i*6+8<<" "<<i*6+9<<" "<<i*6+3<<" "<<i*6+4<<" "<<i*6+10<<" "<<i*6+11<<" "<<i*6+5<<") (1 "<<nDiv<<" "<<ceil((lz-d/2)/delta)<<") simpleGrading (1 1 1)\n";
	}
	ofs << ");\n";
	ofs << "edges (\n";
	ofs << "spline 0 3 (";
	for (int i=0; i < nDiv+1; i++) {
		ofs << "("<<nodes[(nDiv+1)*i][0]<<" "<<nodes[(nDiv+1)*i][1]<<" "<<nodes[(nDiv+1)*i][2]<<") ";
	}
	ofs << ")\n";
	ofs << "spline 4 7 (";
	for (int i=0; i < nDiv+1; i++) {
		ofs << "("<<nodes[(nDiv+1)*i][0]<<" "<<nodes[(nDiv+1)*i][1]<<" "<<nodes[(nDiv+1)*i][2]+d/2<<") ";
	}
	ofs << ")\n";
	ofs << "spline 1 2 (";
	for (int i=0; i < nDiv+1; i++) {
		ofs << "("<<nodes[(nDiv+1)*i+1][0]<<" "<<nodes[(nDiv+1)*i+1][1]<<" "<<nodes[(nDiv+1)*i+1][2]<<") ";
	}
	ofs << ")\n";
	ofs << "spline 5 6 (";
	for (int i=0; i < nDiv+1; i++) {
		ofs << "("<<nodes[(nDiv+1)*i+1][0]<<" "<<nodes[(nDiv+1)*i+1][1]<<" "<<nodes[(nDiv+1)*i+1][2]+d/2<<") ";
	}
	ofs << ")\n";
	for (int j = 2; j < nDiv+1; j++) {
		ofs << "spline "<<6*j<<" "<<6*j+1<<" (";
		for (int i=0; i < nDiv+1; i++) {
			ofs << "("<<nodes[(nDiv+1)*i+j][0]<<" "<<nodes[(nDiv+1)*i+j][1]<<" "<<nodes[(nDiv+1)*i+j][2]<<") ";
		}
		ofs << ")\n";
		ofs << "spline "<<6*j+2<<" "<<6*j+3<<" (";
		for (int i=0; i < nDiv+1; i++) {
			ofs << "("<<nodes[(nDiv+1)*i+j][0]<<" "<<nodes[(nDiv+1)*i+j][1]<<" "<<nodes[(nDiv+1)*i+j][2]+d/2<<") ";
		}
		ofs << ")\n";
	}
	ofs << ");\n";
	ofs << "boundary (\n";
	ofs << "inlet{type patch; faces (";
	ofs << "(3 2 6 7) (7 6 10 11) (2 13 15 6) (6 15 17 10) ";
	for (int i = 2; i < nDiv; i++) {
		ofs << "("<<6*i+1<<" "<<6*i+7<<" "<<6*i+9<<" "<<6*i+3<<") ("<<6*i+3<<" "<<6*i+9<<" "<<6*i+11<<" "<<6*i+5<<") ";
	}
	ofs << ");}\n";
	ofs << "outlet{type patch; faces (";
	ofs << "(0 1 5 4) (4 5 9 8) (1 12 14 5) (5 14 16 9) ";
	for (int i = 2; i < nDiv; i++) {
		ofs << "("<<6*i<<" "<<6*i+6<<" "<<6*i+8<<" "<<6*i+2<<") ("<<6*i+2<<" "<<6*i+8<<" "<<6*i+10<<" "<<6*i+4<<") ";
	}
	ofs << ");}\n";
	ofs << "terrain{type wall; faces (";
	ofs << "(0 1 2 3) (1 12 13 2) ";
	for (int i = 2; i < nDiv; i++) {
		ofs << "("<<6*i<<" "<<6*i+6<<" "<<6*i+7<<" "<<6*i+1<<") ";
	}
	ofs << ");}\n";
	ofs << "top{type patch; faces (";
	ofs << "(8 9 10 11) (9 16 17 10) ";
	for (int i = 2; i < nDiv; i++) {
		ofs << "("<<6*i+4<<" "<<6*i+10<<" "<<6*i+11<<" "<<6*i+5<<") ";
	}
	ofs << ");}\n";
	ofs << "sideE{type patch; faces (("<<6*(nDiv-1)+6<<" "<<6*(nDiv-1)+8<<" "<<6*(nDiv-1)+9<<" "<<6*(nDiv-1)+7<<") ("<<6*(nDiv-1)+8<<" "<<6*(nDiv-1)+10<<" "<<6*(nDiv-1)+11<<" "<<6*(nDiv-1)+9<<"));}\n";
	ofs << "sideW{type patch; faces ((0 4 7 3) (4 8 11 7));});\n";
	ofs << "mergePatchPairs ();";
	ofs.close();
}

void Geometry::threedLayoutRefineComplexTerrain(double d, double h, std::string directory) {
	std::vector< std::vector<double> > turbines;
	std::vector<double> temp(3);
	std::string line;
	std::istringstream iss;
	std::string directoryS = directory+"/turbinesRotated.txt";
	std::ifstream inputs (directoryS.c_str());
	std::getline(inputs, line);
	int rows = 0;
	while (std::getline(inputs, line))
	{
		iss.str(line);
		iss >> temp[0];
		iss >> temp[1];
		iss >> temp[2];
		turbines.push_back(temp);
		rows = rows+1;
		iss.clear();
	}
	inputs.close();
	std::vector<double> slope(rows);
	directoryS = directory+"/terrainRotatedStreamWiseSlope.txt";
	inputs.open(directoryS.c_str());
	std::getline(inputs, line);
	for (int i=0; i<rows; i++)
	{
		std::getline(inputs, line);
		iss.str(line);
		iss >> slope[i];
		iss.clear();
	}
	inputs.close();
	directoryS = directory+"/system/snappyHexMeshDict";
	std::ofstream ofs;
	ofs.open (directoryS.c_str());
	ofs << "FoamFile{version 2.0; format ascii; class dictionary; location \"system\"; object snappyHexMeshDict;}\n";
	ofs << "castellatedMesh true; snap false; addLayers false;\n";
	ofs << "geometry{\n";
	for (int i=0; i<rows; i++)
	{
		ofs << "innerCylinder"<<i+1<<"{type searchableCylinder; ";
		ofs << "point1 ("<<turbines[i][0]<<" "<<turbines[i][1]-d/5<<" "<<turbines[i][2]-slope[i]*d/5+h<<"); point2 ("<<turbines[i][0]<<" "<<turbines[i][1]+d/5<<" "<<turbines[i][2]+slope[i]*d/5+h<<"); radius "<<d/2*1.1<<";}\n";
		ofs << "outerCylinder"<<i+1<<"{type searchableCylinder; ";
		ofs << "point1 ("<<turbines[i][0]<<" "<<turbines[i][1]-d*1.5<<" "<<turbines[i][2]-slope[i]*d*1.5+h<<"); point2 ("<<turbines[i][0]<<" "<<turbines[i][1]+d*1.5<<" "<<turbines[i][2]+slope[i]*d*1.5+h<<"); radius "<<d/2*1.3<<";}\n";
	}
	ofs << "}\n";
	ofs << "castellatedMeshControls{maxLocalCells 10000000; maxGlobalCells 1000000000; minRefinementCells 0; maxLoadUnbalance 0.10; nCellsBetweenLevels 3; features ();\n";
	ofs << "refinementSurfaces{}\n";
	ofs << "resolveFeatureAngle 30; refinementRegions{\n";
	for (int i=0; i<rows; i++)
	{
		ofs << "innerCylinder"<<i+1<<"{mode inside; levels ((2 2));}\n";
		ofs << "outerCylinder"<<i+1<<"{mode inside; levels ((1 1));}\n";
	}
	ofs << "}\n";
	ofs << "locationInMesh ("<<turbines[0][0]<<" "<<turbines[0][1]<<" "<<turbines[0][2]+h<<"); allowFreeStandingZoneFaces true;}\n";
	ofs << "snapControls{nSmoothPatch 4; tolerance 4; nSolveIter 30; nRelaxIter 10;}\n";
	ofs << "meshQualityControls{maxNonOrtho 65; maxBoundarySkewness 20; maxInternalSkewness 4; maxConcave 80; minVol 1e-13; minTetQuality 1e-15; minArea -1; minTwist 0.02;\n";
	ofs << "minDeterminant 0.001; minFaceWeight 0.05; minVolRatio 0.01; minTriangleTwist -1; nSmoothScale 4; errorReduction 0.75; relaxed{maxNonOrtho 75;}}\n";
	ofs << "addLayersControls{relativeSizes true; layers{} expansionRatio 1.0; finalLayerThickness 0.5; minThickness 0.25; nGrow 0; featureAngle 60; nRelaxIter 5;\n";
	ofs << "nSmoothSurfaceNormals 1; nSmoothNormals 3; nSmoothThickness 10; maxFaceThicknessRatio 0.5; maxThicknessToMedialRatio 0.3; minMedianAxisAngle 90;\n";
	ofs << "nBufferCellsNoExtrude 0; nLayerIter 50;}; mergeTolerance 1e-6;\n";
	ofs.close();
}
