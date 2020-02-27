import numpy as np
import sys
from scipy.interpolate import griddata
from operator import add

# Read input files

terrainMeshXList = list()
terrainMeshYList = list()
terrainMeshZList = list()
f = open(sys.argv[1] + "terrainMeshRotated.txt", "r")
f.readline()
for line in f:
	lineParts = line.split()
	terrainMeshXList.append((float(lineParts[0])))
	terrainMeshYList.append((float(lineParts[1])))
	terrainMeshZList.append((float(lineParts[2])))
terrainMeshX = np.array(terrainMeshXList)
terrainMeshY = np.array(terrainMeshYList)
terrainMeshZ = np.array(terrainMeshZList)

# Calculate terrain elevation at turbine positions

turbineXList = list()
turbineYList = list()
f = open(sys.argv[1] + "turbinesRotated.txt", "r")
f.readline()
for line in f:
	lineParts = line.split()
	turbineXList.append((float(lineParts[0])))
	turbineYList.append((float(lineParts[1])))
turbineX = np.array(turbineXList)
turbineY = np.array(turbineYList)

turbineZ = griddata((terrainMeshX,terrainMeshY),terrainMeshZ,(turbineX,turbineY),method="cubic")

# Write new turbine cordinates

f = open(sys.argv[1] + "turbinesRotated.txt", "w")
f.write("X\tY\tZ")
for i in range(len(turbineX)):
	f.write("\n")
	f.write("{:.1f}".format(turbineX[i]))
	f.write("\t")
	f.write("{:.1f}".format(turbineY[i]))
	f.write("\t")
	f.write("{:.3f}".format(turbineZ[i]))

# Calculate terrain slope at turbine positions

delta = 20
turbineXplus = turbineX + delta
turbineYplus = turbineY + delta
turbineXminus = turbineX - delta
turbineYminus = turbineY - delta

turbineZplusX = griddata((terrainMeshX,terrainMeshY),terrainMeshZ,(turbineXplus,turbineY),method="cubic")
turbineZminusX = griddata((terrainMeshX,terrainMeshY),terrainMeshZ,(turbineXminus,turbineY),method="cubic")
turbineZplusY = griddata((terrainMeshX,terrainMeshY),terrainMeshZ,(turbineX,turbineYplus),method="cubic")
turbineZminusY = griddata((terrainMeshX,terrainMeshY),terrainMeshZ,(turbineX,turbineYminus),method="cubic")

terrainStreamWiseSlope = (turbineZplusY-turbineZminusY)/(2*delta)
terrainCrossWiseSlope = (turbineZplusX-turbineZminusX)/(2*delta)

# Write terrains slopes

f = open(sys.argv[1] + "terrainRotatedStreamWiseSlope.txt", "w")
f.write("dZ/dY")
for i in range(len(terrainStreamWiseSlope)):
	f.write("\n")
	f.write("{:.8f}".format(terrainStreamWiseSlope[i]))

f = open(sys.argv[1] + "terrainRotatedCrossWiseSlope.txt", "w")
f.write("dZ/dX")
for i in range(len(terrainCrossWiseSlope)):
	f.write("\n")
	f.write("{:.8f}".format(terrainCrossWiseSlope[i]))
