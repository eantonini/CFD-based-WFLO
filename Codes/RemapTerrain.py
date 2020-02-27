import numpy as np
import sys
from scipy.interpolate import griddata

# Read input files

f = open("InputsPowerCalculation.txt", "r")
for x in range(5):
	f.readline()
line = f.readline()
lineParts = line.split()
temp = float(lineParts[0])
minB = -temp/2
maxB = temp/2
f.readline()
line = f.readline()
lineParts = line.split()
temp = float(lineParts[0])
nDiv = (maxB-minB)/temp

# Extrapolate terrain elevation to the new grid

terrainXList = list()
terrainYList = list()
terrainZList = list()
f = open(sys.argv[1] + "terrainRotated.txt", "r")
f.readline()
for line in f:
	lineParts = line.split()
	terrainXList.append((float(lineParts[0])))
	terrainYList.append((float(lineParts[1])))
	terrainZList.append((float(lineParts[2])))
terrainX = np.array(terrainXList)
terrainY = np.array(terrainYList)
terrainZ = np.array(terrainZList)

n = int(nDiv+1)
xi = yi = np.linspace(minB,maxB,n)
xi,yi = np.meshgrid(xi,yi)
zi = griddata((terrainX,terrainY),terrainZ,(xi,yi),method="linear")

# Smooth the terrain

h = 5 # half width of filter
z_temp = np.zeros((len(zi)+2*h,len(zi[0])+2*h))
z_temp[0:h,0:h] = zi[0][0]*np.ones((h,h))
z_temp[0:h,n+h:n+h+h] = zi[0][n-1]*np.ones((h,h))
z_temp[n+h:n+h+h,0:h] = zi[n-1][0]*np.ones((h,h))
z_temp[n+h:n+h+h,n+h:n+h+h] = zi[n-1][n-1]*np.ones((h,h))
z_temp[0:h,h:n+h] = np.array([zi[0,0:n],]*h)
z_temp[n+h:n+h+h,h:n+h] = np.array([zi[n-1,0:n],]*h)
z_temp[h:n+h,0:h] = np.array([zi[0:n,0],]*h).transpose()
z_temp[h:n+h,n+h:n+h+h] = np.array([zi[0:n,n-1],]*h).transpose()
z_temp[h:n+h,h:n+h] = zi
kernel = 1/(1+2*h)**2*np.ones((1+2*h,1+2*h))
z_smooth = np.zeros((len(zi),len(zi[0])))
for i in range(len(zi)):
    for j in range(len(zi[0])):
        jj = j+h
        ii = i+h
        z_smooth[i][j] = np.sum(z_temp[ii-h:ii+h+1,jj-h:jj+h+1]*kernel)

# Write new terrain coordinates

f = open(sys.argv[1] + "terrainMeshRotated.txt", "w")
f.write("X\tY\tZ")
for i in range(len(xi)):
	for j in range(len(xi[0])):
		f.write("\n")
		f.write("{:.1f}".format(xi[i][j]))
		f.write("\t")
		f.write("{:.1f}".format(yi[i][j]))
		f.write("\t")
		f.write("{:.3f}".format(z_smooth[i][j]))
